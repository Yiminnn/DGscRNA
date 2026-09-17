"""Small, synthetic SLURM probes of the unchanged Python package's data flow."""
from pathlib import Path
from unittest.mock import patch
from contextlib import redirect_stdout, redirect_stderr
import importlib.metadata
import io
import json
import os
import sys
import traceback

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
WORKTREE = ROOT / '.worktrees/reproducibility-gbm-ptc-20260917'
OUT = ROOT / 'results/hvg_ptc_20260916_v1/python_R_audit_20260917'


def run():
    assert os.environ.get('SLURM_JOB_ID'), 'Run scientific probes via SLURM'
    import numpy as np
    import pandas as pd
    import anndata as ad
    import scanpy as sc
    import torch
    sys.path.insert(0, str(WORKTREE))
    from dgscrna.core import preprocessing, clustering, utils, deep_learning
    torch.set_num_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    reports = []

    def record(name, fn):
        log = io.StringIO()
        try:
            with redirect_stdout(log), redirect_stderr(log):
                result = fn()
            reports.append(dict(probe=name, executed=True, result=result))
        except Exception:
            reports.append(dict(probe=name, executed=False, error=traceback.format_exc()))
        (OUT / (name + '.txt')).write_text(log.getvalue())

    def scaled_deg():
        rng = np.random.default_rng(7)
        x = rng.poisson(2, (120, 600)).astype(np.float64)
        x[:60, :20] += 25
        x[60:, 20:40] += 25
        a = ad.AnnData(x)
        a.var_names = [f'G{i}' for i in range(600)]
        a.obs['truth_block'] = pd.Categorical(['A']*60 + ['B']*60)
        a = preprocessing.preprocess_adata(a, min_genes=1, n_pcs=10, n_neighbors=10)
        clustering.find_markers(a, groupby='truth_block', n_genes=100)
        lfc = a.uns['rank_genes_groups']['logfoldchanges']
        values = np.concatenate([lfc[c] for c in lfc.dtype.names])
        return dict(lognorm_layer_saved='lognorm' in a.layers, mitochondrial_metric_computed='pct_counts_mt' in a.obs,
                    n_hvg=int(a.var.highly_variable.sum()), n_logfc=len(values), n_nonfinite_logfc=int((~np.isfinite(values)).sum()),
                    failure='preprocess scales X without preserving lognorm; default DEG falls back to scaled X')

    def multiple_partitions():
        a = ad.AnnData(np.ones((12, 4)))
        a.obsm['X_pca'] = np.ones((12, 2))
        captured = []

        def partition(a, methods, **kwargs):
            for method in methods:
                a.obs[method+'_clusters'] = [method+'_A']*6 + [method+'_B']*6
            return a

        def markers(a, groupby, **kwargs):
            a.uns['rank_genes_groups'] = {'originating_partition':groupby}
            return a

        def score(a, marker_sets, cluster_key, **kwargs):
            captured.append(dict(scored_partition=cluster_key,
                                 deg_partition=a.uns['rank_genes_groups']['originating_partition']))
            return a

        with patch.object(utils, 'run_clustering', partition), patch.object(utils, 'find_markers', markers), \
             patch.object(utils, 'load_marker_sets', return_value={'panel': {'type':['G0']}}), \
             patch.object(utils, 'score_cell_types', score):
            utils.run_dgscrna_pipeline(a, 'synthetic', clustering_methods=['leiden','hdbscan'], use_deep_learning=False)
        return dict(calls=captured, all_partitions_use_own_DEG=all(r['scored_partition']==r['deg_partition'] for r in captured),
                    probe_type='mocked data-flow isolation; no clustering algorithm performance measured')

    def harmony_routing():
        parts = []
        for group in range(2):
            a = ad.AnnData(np.ones((6, 4)))
            a.obs_names = [f'{group}_{i}' for i in range(6)]
            a.obsm['X_pca'] = np.arange(12, dtype=float).reshape(6, 2) + group
            parts.append(a)

        def harmony(a, key, **kwargs):
            a.obsm['X_pca_harmony'] = a.obsm['X_pca'] + 100

        captured = {}

        class CaptureHDBSCAN:
            def __init__(self, **kwargs):
                pass
            def fit_predict(self, x):
                captured['matrix'] = x.copy()
                return np.zeros(len(x), dtype=int)

        with patch.object(sc.external.pp, 'harmony_integrate', harmony):
            a = preprocessing.integrate_datasets(parts, method='harmony')
        a.uns['neighbors'] = {'probe_only':True}
        with patch.object(clustering.hdbscan, 'HDBSCAN', CaptureHDBSCAN):
            clustering.run_clustering(a, methods=['hdbscan'], cluster_space='pca')
        return dict(corrected_basis_exists='X_pca_harmony' in a.obsm,
                    hdbscan_consumes_corrected_basis=bool(np.array_equal(captured['matrix'], a.obsm['X_pca_harmony'])),
                    hdbscan_consumes_uncorrected_basis=bool(np.array_equal(captured['matrix'], a.obsm['X_pca'])),
                    probe_type='mocked Harmony and clustering; verifies representation routing only')

    def unknown_training():
        a = ad.AnnData(np.ones((60, 8)))
        a.obs['initial'] = ['TypeA']*20 + ['Unknown']*20 + ['Undecided']*20
        *_, mapping = deep_learning.prepare_training_data(a, 'initial', use_highly_variable=False)
        return dict(training_classes=list(mapping.values()), Unknown_is_a_training_class='Unknown' in mapping.values(),
                    Undecided_is_a_training_class='Undecided' in mapping.values())

    record('default_preprocessing_DEG', scaled_deg)
    record('multiple_partition_DEG_isolation', multiple_partitions)
    record('harmony_corrected_basis_routing', harmony_routing)
    record('Unknown_training_class', unknown_training)
    versions = {}
    for package in ['numpy','scipy','pandas','scanpy','anndata','torch','hdbscan','harmonypy','scikit-learn','umap-learn']:
        versions[package] = importlib.metadata.version(package)
    report = dict(job=os.environ['SLURM_JOB_ID'], package_commit='4bf17c4cb9518ab427fa92e8bf5bb4d594d30e00',
                  package_path=str(WORKTREE / 'dgscrna'), inputs='synthetic; no patient data or production outputs loaded',
                  package_modified=False, versions=versions, probes=reports)
    (OUT / 'package_probe_results.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)
    assert all(r['executed'] for r in reports), 'An audit probe failed to execute; inspect report'


if __name__ == '__main__':
    run()
