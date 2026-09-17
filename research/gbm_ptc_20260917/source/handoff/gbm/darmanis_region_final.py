"""Darmanis donor×Location独立终端注释；固定参数，不救援小region。

复用已核验的panel读取、strict7映射、counts/QC与评分规则。
DL只补非Noise未知细胞；全filtered scaled基因是唯一实际特征矩阵。
"""
import argparse
from datetime import datetime, timezone
import importlib.metadata
import inspect
import json
import os
from pathlib import Path
import traceback

import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_fscore_support

from darmanis_validate_20260910 import (
    ABSTAIN, GOLD, MARKERS, ROOT, SETS, VOCAB, digest, evaluate,
    load_panels, map_call, qc_mask,
)

UNITS = [('BT_S4', 'Tumor'), ('BT_S4', 'Periphery'),
         ('BT_S2', 'Tumor'), ('BT_S2', 'Periphery'),
         ('BT_S1', 'Tumor'), ('BT_S1', 'Periphery'),
         ('BT_S6', 'Tumor'), ('BT_S6', 'Periphery'), ('BT_S6', 'Distant')]
OUTROOT = ROOT / 'results/darmanis_region_v1/annotation_v1'


def select_unit(a, donor, location):
    mask = (a.obs['donor_id'].astype(str) == donor) & (a.obs['Location'].astype(str) == location)
    result = a[mask.values].copy()
    if not result.n_obs or not result.obs_names.is_unique:
        raise ValueError(f'不存在的unit或重复cell ID: {donor}/{location}')
    return result


def oracle_bounds(gold, reachable):
    present = set(gold) & set(VOCAB)
    reachable = set(reachable) & set(VOCAB)
    return {'oracle_present_upper_bound': len(present & reachable) / len(present),
            'oracle_fixed7_sample_upper_bound': len(present & reachable) / len(VOCAB),
            'oracle_fixed7_vocab_upper_bound': len(reachable) / len(VOCAB)}


def refine_final(a, annotation_key, train_fn=None, predict_fn=None):
    from dgscrna.core.deep_learning import train_deep_model, predict_cell_types
    train_fn = train_fn or train_deep_model
    predict_fn = predict_fn or predict_cell_types
    raw = a.obs[annotation_key].astype(str)
    noise = a.obs['cluster'].astype(str).eq('Noise')
    pool = ~noise & raw.isin(ABSTAIN)
    training = ~noise & ~raw.isin(ABSTAIN)
    classes = sorted(raw[training].unique().tolist())
    final = raw.copy().astype(object)
    final.loc[noise | pool] = 'Unknown'
    lineage = pd.Series('marker_retained', index=a.obs_names, dtype=object)
    lineage.loc[noise] = 'noise_terminal'
    lineage.loc[pool] = 'unresolved_no_training'
    status = {'n_nonnoise': int((~noise).sum()), 'n_noise': int(noise.sum()),
              'n_pool': int(pool.sum()), 'n_training': int(training.sum()),
              'n_training_classes': len(classes), 'training_classes': classes,
              'training_attempted': False, 'training_executed': False,
              'prediction_executed': False, 'final_valid': True}
    if not pool.any():
        status['dl_status'] = 'empty_pool'
    elif len(classes) < 2:
        status.update(dl_status='lt2_training_classes', final_valid=False)
    else:
        subset = a[~noise.values].copy()
        labels = raw.loc[subset.obs_names].astype(object)
        labels.loc[labels.isin(ABSTAIN)] = 'Undecided'
        subset.obs['_dl_seed_panel'] = labels
        assert subset.obsm['X_scaled'].shape == (subset.n_obs, a.n_vars)
        try:
            status['training_attempted'] = True
            model = train_fn(subset, '_dl_seed_panel', epochs=15, batch_size=256,
                             use_highly_variable=False, device='cpu', random_state=42)
            status['training_executed'] = True
            if model['num_classes'] != len(classes) or model['input_dim'] != a.n_vars:
                raise ValueError('DL实际类别数或features维数与manifest不符')
            result = predict_fn(model, subset, '_dl_seed_panel', probability_threshold=.9,
                                use_highly_variable=False)
            if not result.index.equals(subset.obs_names):
                raise ValueError('DL返回cell ID错序或缺失')
            if not set(result.astype(str)) <= set(classes) | ABSTAIN:
                raise ValueError('DL返回未知panel')
            if not np.array_equal(result.loc[training[training].index].astype(str), raw.loc[training]):
                raise ValueError('DL修改了已确定的panel')
            result = result.astype(object)
            result.loc[result.isin(ABSTAIN)] = 'Unknown'
            final.loc[subset.obs_names] = result
            lineage.loc[pool] = np.where(final.loc[pool].eq('Unknown'),
                                         'dl_low_confidence', 'dl_confident')
            status.update(dl_status='dl_completed', prediction_executed=True)
        except Exception as exc:
            status.update(dl_status='dl_error', final_valid=False,
                          error_type=type(exc).__name__, error=str(exc),
                          traceback=traceback.format_exc())
    status['n_dl_assigned'] = int(lineage.eq('dl_confident').sum())
    status['lineage_counts'] = lineage.value_counts().to_dict()
    return final, lineage, status


def save_manifest(out, manifest):
    tmp = out / 'manifest.json.tmp'
    tmp.write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    tmp.replace(out / 'manifest.json')


def compute(out, donor, location, manifest):
    import anndata as ad
    import hdbscan
    import scanpy as sc
    import scipy.sparse as sp
    import torch
    from dgscrna.core.deep_learning import train_deep_model, predict_cell_types
    from dgscrna.core.marker_scoring import score_cell_types

    torch.set_num_threads(1)
    src = ROOT / 'data_bench/brain_GBM/brain_GBM.h5ad'
    sources = [src, Path(__file__), Path(__file__).with_name('test_darmanis_region_final.py'),
               Path(__file__).with_name('darmanis_region_final.sbatch'),
               Path(__file__).with_name('darmanis_validate_20260910.py'),
               Path(inspect.getfile(score_cell_types)), Path(inspect.getfile(train_deep_model)),
               ROOT/'DGscRNA/dgscrna/models/deep_model.py', MARKERS/'mapping_L1_v3.csv']
    sources += [MARKERS / f'{name}.csv' for name in SETS]
    manifest.update(
        dependency_versions={name: importlib.metadata.version(name) for name in
                             ['dgscrna', 'scanpy', 'hdbscan', 'anndata', 'numpy', 'scipy',
                              'scikit-learn', 'umap-learn', 'torch']},
        input_source=str(src), source_hashes={str(p): digest(p) for p in sources},
        params={'marker_sets': list(SETS), 'cutoff': 'none', 'seed': 42,
                'hdbscan_min_cluster_size': 15, 'hdbscan_min_samples': 15,
                'gene_min_cells': 3, 'normalize_total': 1e4, 'scale_max_value': 10,
                'pca_n_comps': 30, 'neighbors': 15, 'umap_n_components': 2,
                'deg_method': 'wilcoxon', 'deg_top_n': 100, 'deg_logfc_gt': 1.,
                'deg_padj_lt': .05, 'deg_layer': 'lognorm', 'deg_reference': 'rest_including_noise',
                'dl_epochs': 15, 'dl_batch_size': 256, 'dl_probability_threshold': .9,
                'dl_use_highly_variable': False, 'dl_device': 'cpu'},
        qc_policy='author_qc_plus_nGene200_visibleMT15',
        qc_note='沿用发布的作者QC细胞；额外nGene>=200，有可见MT时要求MT<=15%；不额外DoubletFinder。',
        feature_policy='全部region内min_cells>=3的基因；共同scale X显式写入X_scaled；无HVG筛选。',
        counts_provenance='CELLxGENE下载X经全值counts-like验证；不是原始全部测序细胞。',
        reference_note='作者cell_type经预设strict7映射，仅称concordance，不是独立实测accuracy。',
        oracle_note='fixed7_sample为可达∩本unit present / 7；present为可达∩present / present数。',
        stage='load_counts')
    save_manifest(out, manifest)
    sets, mapping = load_panels()
    a = select_unit(ad.read_h5ad(src), donor, location)
    values = a.X.data if sp.issparse(a.X) else np.asarray(a.X).ravel()
    if not (np.isfinite(values).all() and (values >= 0).all() and
            np.isclose(values, np.rint(values), atol=1e-6, rtol=0).all()):
        raise ValueError('X不是非负有限counts-like矩阵，禁止继续归一化')
    a.var['source_id'] = a.var_names.astype(str)
    symbols = a.var['feature_name'].astype(str)
    if symbols.duplicated().any() or symbols.isin(['', 'nan']).any():
        raise ValueError('重复/缺失gene symbol')
    a.var_names = symbols.values
    a.var['mt'] = a.var_names.str.upper().str.startswith('MT-')
    has_mt = bool(a.var['mt'].any())
    sc.pp.calculate_qc_metrics(a, qc_vars=['mt'] if has_mt else [], percent_top=None,
                              log1p=False, inplace=True)
    if not has_mt:
        a.obs['pct_counts_mt'] = np.nan
    qc = a.obs[['donor_id', 'Location', 'Selection', 'cell_type', 'total_counts',
                'n_genes_by_counts', 'pct_counts_mt']].copy()
    qc.index.name = 'cell'
    qc['nGene200_pass'] = qc['n_genes_by_counts'] >= 200
    qc['visible_MT_available'] = has_mt
    qc['visible_MT15_pass'] = (qc['pct_counts_mt'] <= 15) if has_mt else pd.NA
    qc['retained'] = qc_mask(qc, has_mt, 'author_qc')
    qc['filtered_by_current_qc'] = ~qc['retained']
    qc.to_csv(out / 'qc_cells.csv')
    manifest.update(input_shape=list(a.shape), n_input=a.n_obs,
                    n_retained=int(qc.retained.sum()), n_filtered=int((~qc.retained).sum()),
                    mt_qc='computed' if has_mt else 'unavailable_in_released_expression_matrix',
                    n_visible_mt_genes=int(a.var['mt'].sum()), stage='preprocess')
    a = a[qc.retained.values].copy()
    if a.n_obs <= 31:
        raise ValueError('QC后细胞不足固定PCA30；不调参救援')
    sc.pp.filter_genes(a, min_cells=3)
    if a.n_vars <= 30:
        raise ValueError('filtered基因不足固定PCA30；不调参救援')
    a.var[['source_id']].rename_axis('symbol').to_csv(out / 'genes_filtered.csv')
    manifest.update(qc_shape=list(a.shape), actual_dl_n_features=a.n_vars)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.layers['lognorm'] = a.X.copy()
    sc.pp.scale(a, max_value=10)
    a.obsm['X_scaled'] = np.asarray(a.X, dtype=np.float32).copy()
    sc.tl.pca(a, n_comps=30, svd_solver='arpack', random_state=42, mask_var=None)
    sc.pp.neighbors(a, n_neighbors=15, use_rep='X_pca', random_state=42)
    sc.tl.umap(a, n_components=2, random_state=42)
    label = hdbscan.HDBSCAN(min_cluster_size=15, min_samples=15).fit_predict(a.obsm['X_umap'])
    cl = np.asarray(['Noise' if i == -1 else str(i) for i in label])
    a.obs['cluster'] = pd.Categorical(cl)
    groups = sorted(set(cl) - {'Noise'})
    np.savez_compressed(out / 'embedding.npz', pca=a.obsm['X_pca'], umap=a.obsm['X_umap'],
                        cells=np.asarray(a.obs_names, dtype=str), cluster=cl)
    gold = a.obs['cell_type'].map(GOLD)
    if gold.isna().any():
        raise ValueError('未定义的作者cell_type')
    pd.DataFrame({'cell': a.obs_names, 'donor': donor, 'Location': location,
                  'gold': gold.values, 'cluster': cl}).to_csv(out / 'cells_final.csv', index=False)
    deg_available = len(groups) >= 2
    manifest.update(stage='deg', n_clusters=len(groups), n_noise=int((cl == 'Noise').sum()),
                    cluster_counts=pd.Series(cl).value_counts().to_dict(),
                    deg_status='computed' if deg_available else 'insufficient_clusters')
    save_manifest(out, manifest)
    if deg_available:
        sc.tl.rank_genes_groups(a, 'cluster', groups=groups, reference='rest', method='wilcoxon',
                               n_genes=min(100, a.n_vars), layer='lognorm', use_raw=False)
        deg = sc.get.rank_genes_groups_df(a, group=None)
        deg.to_csv(out / 'deg.csv', index=False)
        significant = deg[(deg.pvals_adj < .05) & (deg.logfoldchanges > 1.)]
        manifest['significant_deg_counts'] = significant.groupby('group', observed=True).size().to_dict()
    coverage, rows, perclass, predictions = [], [], [], []
    manifest['marker_status'] = {}
    for name, panels in sets.items():
        reachable = {map_call(name, p, mapping) for p in panels} & set(VOCAB)
        for panel, genes in panels.items():
            coverage.append({'marker_set': name, 'panel': panel, 'n_panel_genes': len(genes),
                             'n_detected_unique': len(set(genes) & set(a.var_names)),
                             'target': map_call(name, panel, mapping)})
        manifest.update(stage='marker_and_dl', current_marker=name)
        save_manifest(out, manifest)
        key = f'{name}_cluster_none'
        if deg_available:
            score_cell_types(a, panels, 'cluster', marker_set_name=name, cutoffs=['none'],
                             min_logfc=1., max_pval_adj=.05, random_state=42)
        else:
            a.obs[key] = 'Unknown'
        final, lineage, state = refine_final(a, key)
        if not deg_available:
            state.update(dl_status='all_noise' if not groups else 'insufficient_clusters',
                         final_valid=not bool(groups))
            lineage.loc[cl != 'Noise'] = 'unresolved_deg_unavailable'
            state['lineage_counts'] = lineage.value_counts().to_dict()
        pred = np.asarray([map_call(name, p, mapping) for p in final])
        scalar = {k: v for k, v in state.items() if not isinstance(v, (dict, list)) and k != 'traceback'}
        common = {'donor': donor, 'Location': location, 'unit': f'{donor}__{location}',
                  'marker_set': name, 'cutoff': 'none', 'stage': 'final', 'n_clusters': len(groups)}
        metric = {**common, **evaluate(gold.values, pred), **oracle_bounds(gold, reachable), **scalar}
        rows.append(metric)
        precision, recall, f1, support = precision_recall_fscore_support(
            gold.values, pred, labels=list(VOCAB), zero_division=0)
        for i, cellclass in enumerate(VOCAB):
            perclass.append({**common, 'class': cellclass, 'precision': precision[i],
                             'recall': recall[i], 'f1': f1[i], 'support': int(support[i]),
                             'present': bool(support[i]), 'final_valid': state['final_valid']})
        predictions.append(pd.DataFrame({**common, 'cell': a.obs_names, 'gold': gold.values,
                                        'cluster': cl, 'final_panel': final.values,
                                        'final_label': pred, 'prediction_lineage': lineage.values,
                                        'dl_status': state['dl_status'], 'final_valid': state['final_valid']}))
        manifest['marker_status'][name] = state
        print(json.dumps({'marker_set': name, **state}, ensure_ascii=False), flush=True)
        save_manifest(out, manifest)
    pd.concat(predictions, ignore_index=True).to_csv(out / 'final_predictions.csv', index=False)
    pd.DataFrame(rows).to_csv(out / 'metrics_final.csv', index=False)
    pd.DataFrame(perclass).to_csv(out / 'perclass_final.csv', index=False)
    pd.DataFrame(coverage).to_csv(out / 'panel_coverage.csv', index=False)
    if any(not state['final_valid'] for state in manifest['marker_status'].values()):
        manifest['status'] = 'computed_with_fallback_not_reviewed'
    else:
        manifest['status'] = 'computed_not_reviewed'
    manifest.update(stage='final', n_output_cells=a.n_obs,
                    finished_at=datetime.now(timezone.utc).isoformat(),
                    output_hashes={p.name: digest(p) for p in sorted(out.iterdir())
                                   if p.is_file() and p.name != 'manifest.json'})
    save_manifest(out, manifest)
    (out / 'COMPLETE').write_text(manifest['status'] + '\n')


def run(args):
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('科学计算仅允许在SLURM作业内运行')
    donor, location = UNITS[args.unit_index]
    out = OUTROOT / f'{donor}__{location}'
    out.mkdir(parents=True, exist_ok=False)
    manifest = {'donor': donor, 'Location': location, 'unit': out.name, 'status': 'started',
                'started_at': datetime.now(timezone.utc).isoformat(),
                'slurm_job': os.environ['SLURM_JOB_ID'],
                'slurm_array_job': os.environ.get('SLURM_ARRAY_JOB_ID'),
                'slurm_array_task': os.environ.get('SLURM_ARRAY_TASK_ID')}
    save_manifest(out, manifest)
    try:
        compute(out, donor, location, manifest)
    except BaseException as exc:
        manifest.update(status='failed', error_type=type(exc).__name__, error=str(exc),
                        traceback=traceback.format_exc())
        save_manifest(out, manifest)
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--unit-index', type=int, choices=range(len(UNITS)), required=True)
    run(parser.parse_args())
