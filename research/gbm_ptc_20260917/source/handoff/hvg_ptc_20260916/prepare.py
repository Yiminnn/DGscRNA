#!/usr/bin/env python3
"""Prepare label-free, frozen-cell inputs once per sample, on SLURM only."""
import argparse
import json
import os
from pathlib import Path
import resource
import sys
import time
import traceback

from common import (ROOT, OUT, INPUTS, FEATURES, EXTRA_FEATURES, MARKERS, MARKER,
                    require_slurm, runtime_record, sha, utc, version_record, write_json)


def prepare(sample):
    require_slurm()
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import sparse
    sys.path.insert(0, str(ROOT / 'handoff/g274_table4'))
    from fit_runner import load_clean_counts

    dest = OUT / 'prepared' / sample
    dest.mkdir(parents=True, exist_ok=True)
    if (dest / 'manifest.json').exists() and not (dest / 'PREPARED').exists():
        previous = json.loads((dest / 'manifest.json').read_text())
        if previous.get('status') == 'failed':
            archive = dest / f'attempt_{previous.get("slurm_job_id", "unknown")}.json'
            if not archive.exists():
                archive.write_text((dest / 'manifest.json').read_text())
    source = INPUTS / sample
    frozen = json.loads((source / 'manifest.json').read_text())
    expected = frozen['artifacts']['counts_gene_filtered.h5ad']['sha256']
    started = time.perf_counter()
    meta = {'sample': sample, 'status': 'running', 'started_at': utc(),
            'source': str(source), 'source_counts_sha256': expected,
            'gold_used_for_preprocessing': False, 'versions': version_record(),
            'code': {str(p): sha(p) for p in [Path(__file__), Path(__file__).with_name('common.py')]},
            **runtime_record()}
    if (dest / 'PREPARED').exists():
        old = json.loads((dest / 'manifest.json').read_text())
        if old['status'] == 'completed' and old['source_counts_sha256'] == expected and old['code'] == meta['code']:
            for name, h in old['outputs'].items():
                if not (dest / name).is_file() or sha(dest / name) != h:
                    raise ValueError(f'Prepared artifact failed checksum: {sample}/{name}')
            print(f'{sample}: verified cached preparation', flush=True)
            return
        raise ValueError('Prepared input differs; preserve previous version and investigate')
    write_json(dest / 'manifest.json', meta)
    try:
        if sha(source / 'counts_gene_filtered.h5ad') != expected:
            raise ValueError('Frozen input checksum mismatch')
        a = load_clean_counts(source)
        counts = sparse.csr_matrix(a.X, dtype=np.float32)
        values = counts.data
        if not (np.isfinite(values).all() and (values >= 0).all() and np.equal(values, np.rint(values)).all()):
            raise ValueError('Input is not nonnegative integer counts')
        (dest / 'cells.tsv').write_text('\n'.join(a.obs_names) + '\n')
        (dest / 'genes.tsv').write_text('\n'.join(a.var_names) + '\n')
        meta.update(n_cells=a.n_obs, n_genes=a.n_vars, stage='HVG', feature_sets={})
        write_json(dest / 'manifest.json', meta)
        selected = {'all': np.arange(a.n_vars, dtype=np.int32)}
        # VST requires counts. No reference labels or annotations enter this object.
        for n in [500, 1000, 2000, 3000, 5000]:
            if a.n_vars < n:
                raise ValueError(f'Requested {n} genes but only {a.n_vars} are available')
            try:
                sc.pp.highly_variable_genes(a, n_top_genes=n, flavor='seurat_v3', subset=False)
            except ValueError as exc:
                if 'singular' not in str(exc).lower():
                    raise
                # Preserve a failed VST condition without changing span or fabricating a gene set.
                # All-gene and independent legacy-dispersion controls remain computable.
                meta.setdefault('unavailable_features', {})[f'hvg{n}'] = {
                    'status': 'structural_failure', 'stage': 'VST_loess',
                    'error': str(exc), 'span': 0.3, 'flavor': 'seurat_v3'}
                continue
            idx = np.flatnonzero(a.var['highly_variable'].to_numpy()).astype(np.int32)
            if len(idx) != n:
                raise ValueError('VST selected a different number of genes')
            selected[f'hvg{n}'] = idx
            cols = [x for x in ['highly_variable', 'highly_variable_rank', 'means',
                                'variances', 'variances_norm'] if x in a.var]
            a.var[cols].rename_axis('gene').to_csv(dest / f'hvg{n}_ranking.csv.gz')
        a.X = counts.copy()
        sc.pp.normalize_total(a, target_sum=10000)
        sc.pp.log1p(a)
        sparse.save_npz(dest / 'lognorm.npz', sparse.csr_matrix(a.X), compressed=True)
        for n in [2000, 5000]:
            sc.pp.highly_variable_genes(a, n_top_genes=n, flavor='seurat', subset=False)
            selected[f'seurat{n}'] = np.flatnonzero(a.var['highly_variable'].to_numpy()).astype(np.int32)
        panel = pd.read_csv(MARKERS / f'{MARKER}.csv', index_col=0)
        marker_genes = set(panel.stack().dropna().astype(str))
        marker_idx = np.flatnonzero(a.var_names.isin(marker_genes))
        if 'hvg2000' in selected:
            selected['hvg2000_markers'] = np.union1d(selected['hvg2000'], marker_idx).astype(np.int32)
        else:
            meta.setdefault('unavailable_features', {})['hvg2000_markers'] = {
                'status': 'structural_failure', 'stage': 'VST_loess',
                'error': 'Requires unavailable hvg2000 mask'}
        coverage = []
        for name, idx in selected.items():
            np.save(dest / f'indices_{name}.npy', idx, allow_pickle=False)
            feature_genes = set(a.var_names[idx])
            meta['feature_sets'][name] = {'n_features': len(idx),
                'selection': 'all min_cells>=3' if name == 'all' else
                    ('Seurat log-dispersion' if name.startswith('seurat') else
                     'VST2000 union frozen marker genes' if name.endswith('_markers') else 'Seurat v3 VST on counts'),
                'indices_sha256': sha(dest / f'indices_{name}.npy')}
            for col in panel:
                genes = set(panel[col].dropna().astype(str))
                detected = genes.intersection(a.var_names)
                coverage.append({'feature': name, 'panel': col, 'n_original': len(genes),
                                 'n_detected': len(detected), 'n_in_geometry': len(genes & feature_genes)})
        pd.DataFrame(coverage).to_csv(dest / 'marker_coverage.csv', index=False)
        meta['stage'] = 'scale_all_genes'
        write_json(dest / 'manifest.json', meta)
        a.X = a.X.toarray()
        sc.pp.scale(a, max_value=10)
        scaled = np.asarray(a.X, dtype=np.float32)
        if not np.isfinite(scaled).all():
            raise ValueError('Scaled expression contains nonfinite values')
        np.save(dest / 'scaled_all.npy', scaled, allow_pickle=False)
        old = ROOT / 'results/g274_cohort/fit_v1' / sample / 'prepare/X.npy'
        meta['old_allgene_scaled_match'] = sha(old) == sha(dest / 'scaled_all.npy') if old.exists() else None
        meta.update(status='completed', stage='complete', finished_at=utc(),
                    elapsed_seconds=time.perf_counter() - started,
                    peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
                    outputs={p.name: sha(p) for p in dest.iterdir()
                             if p.is_file() and p.name not in ['manifest.json', 'PREPARED']
                             and not p.name.startswith('attempt_')})
        write_json(dest / 'manifest.json', meta)
        (dest / 'PREPARED').write_text(sha(dest / 'manifest.json') + '\n')
        print(json.dumps({'sample': sample, 'status': 'completed', 'seconds': meta['elapsed_seconds'],
                          'shape': list(scaled.shape), 'old_match': meta['old_allgene_scaled_match']}), flush=True)
    except BaseException as exc:
        meta.update(status='failed', error_type=type(exc).__name__, error=str(exc),
                    traceback=traceback.format_exc(), finished_at=utc())
        write_json(dest / 'manifest.json', meta)
        raise


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--sample')
    p.add_argument('--sample-list', type=Path)
    p.add_argument('--index', type=int)
    args = p.parse_args()
    sample = args.sample or args.sample_list.read_text().split()[args.index]
    prepare(sample)
