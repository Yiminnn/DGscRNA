"""固定 Darmanis query 的 marker-only 验证；不根据分数换样本或默认参数。

QC策略显式指定：strict_mt要求MT<=15%；author_qc在发布矩阵无MT时记为不可评估，
沿用作者已QC细胞，并加nFeature>=200、gene min_cells=3；不是完整重做QC。
不筛HVG、不删高表达top基因；PCA30 -> kNN15 -> UMAP2 -> HDBSCAN mcs15/ms15。
mcs50 为预先固定的包默认参数敏感性，不取两者最优作为主结果。
DEG：全基因lognorm上的Wilcoxon，保留包对齐top100，p_adj<0.05、logFC>1。
调用包score_cell_types：原panel分母、singleton penalty、ties、mean/0.5/none。
主cutoff为mean。输出marker-only，绝不标为含DL的完整DG-scRNA。
Smart-seq2不套用10x multiplet-rate；本轮不额外用DoubletFinder删细胞。
counts-like X的全值检查是必要条件，不等同于独立核实原作者的全部处理历史。
"""
import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import importlib.metadata
import inspect
import csv
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import f1_score, precision_recall_fscore_support

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
MARKERS = ROOT / 'handoff/markers_v3'
SETS = ('CM2_glioma_other', 'CM2_glioblastoma', 'UNION_CellMarker2_brain', 'UNION_all',
        'BrainAtlas112', 'CARE_TME', 'Liu_devbrain', 'CM2_brain_normal')
VOCAB = ('Malignant', 'Myeloid', 'Astrocyte', 'Oligodendrocyte', 'OPC', 'Neuron', 'Vascular')
GOLD = {'neoplastic cell':'Malignant', 'myeloid cell':'Myeloid', 'astrocyte':'Astrocyte',
        'oligodendrocyte':'Oligodendrocyte', 'oligodendrocyte precursor cell':'OPC',
        'neuron':'Neuron', 'vascular lymphangioblast':'Vascular'}
BRIDGE = {'Malignant':'Malignant', 'TAM':'Myeloid', 'Astrocyte':'Astrocyte',
          'Oligodendrocyte':'Oligodendrocyte', 'OPC':'OPC', 'Endothel':'Vascular',
          'Pericyte':'Vascular', 'Excitatory neuron':'Neuron', 'Inhibitory neuron':'Neuron',
          'AMBIGUOUS_NEURON':'Neuron'}
ABSTAIN = {'Unknown', 'Undecided', 'Noise', 'noise'}
DONORS = ('BT_S4', 'BT_S2', 'BT_S1', 'BT_S6')  # 按输入细胞数降序预设，独立于注释分数。


def load_panels():
    mapping = {}
    with (MARKERS / 'mapping_L1_v3.csv').open() as f:
        for row in csv.DictReader(f):
            key = (row['marker_set'], row['panel'])
            if key in mapping:
                raise ValueError(f'重复映射键: {key}')
            mapping[key] = row['gold_class']
    sets = {}
    for name in SETS:
        frame = pd.read_csv(MARKERS / f'{name}.csv', index_col=0)
        sets[name] = {}
        for panel in frame:
            genes = [str(g).strip() for g in frame[panel].dropna() if str(g).strip()]
            if not genes:
                raise ValueError(f'空panel: {name}/{panel}')
            if (name, panel) not in mapping:
                raise ValueError(f'缺少精确映射: {name}/{panel}')
            sets[name][panel] = genes
    return sets, mapping


def map_call(name, panel, mapping):
    if panel in ABSTAIN:
        return 'Undecided'
    if (name, panel) not in mapping:
        raise ValueError(f'未定义的预测panel: {name}/{panel}')
    return BRIDGE.get(mapping[(name, panel)], 'OFFVOCAB')


def evaluate(gold, prediction):
    gold = np.asarray(gold, dtype=str)
    pred = np.asarray(prediction, dtype=str)
    if gold.size == 0 or gold.shape != pred.shape or not set(gold) <= set(VOCAB):
        raise ValueError('gold/pred长度或gold词表非法')
    ab = np.isin(pred, list(ABSTAIN))
    off = ~np.isin(pred, list(VOCAB)) & ~ab
    present = [c for c in VOCAB if c in gold]
    return {'n_cells': int(gold.size), 'n_present_classes': len(present),
            'macro_f1_fixed7': float(f1_score(gold, pred, labels=list(VOCAB), average='macro', zero_division=0)),
            'macro_f1_present': float(f1_score(gold, pred, labels=present, average='macro', zero_division=0)),
            'concordance': float(np.mean(gold == pred)), 'abstain_rate': float(ab.mean()),
            'offvocab_rate': float(off.mean()), 'coverage': float(1-ab.mean())}


def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


@contextmanager
def record_run(out, manifest):
    out.mkdir(parents=True, exist_ok=False)
    def save():
        temp = out/'manifest.json.tmp'
        temp.write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
        temp.replace(out/'manifest.json')
    manifest.update(status='started', started_at=datetime.now(timezone.utc).isoformat())
    save()
    try:
        yield manifest
    except BaseException as exc:
        manifest.update(status='failed', error_type=type(exc).__name__, error=str(exc))
        save()
        raise
    else:
        manifest.update(status='computed_not_reviewed', finished_at=datetime.now(timezone.utc).isoformat())
        save()
        (out/'COMPLETE').write_text('marker-only computed; review and figures pending\n')


def qc_mask(qc, has_mt, policy):
    if policy not in ('strict_mt', 'author_qc'):
        raise ValueError('未知QC策略')
    if not has_mt and policy == 'strict_mt':
        raise ValueError('未识别MT基因，不能执行严格MT质控')
    mask = qc['n_genes_by_counts'] >= 200
    if has_mt:
        mask &= qc['pct_counts_mt'].notna() & (qc['pct_counts_mt'] <= 15)
    return mask


def run(args):
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('矩阵计算只允许在SLURM作业内运行')
    out = ROOT / 'results/darmanis_validation_20260910' / args.run_id / args.donor
    with record_run(out, {'donor':args.donor, 'qc_policy':args.qc_policy,
                          'slurm_job':os.environ['SLURM_JOB_ID']}) as manifest:
        compute(args, out, manifest)


def compute(args, out, manifest):
    import anndata as ad
    import hdbscan
    import scanpy as sc
    import scipy.sparse as sp
    from dgscrna.core.marker_scoring import score_cell_types

    src = ROOT / 'data_bench/brain_GBM/brain_GBM.h5ad'
    scoring_source = Path(inspect.getfile(score_cell_types))
    files = [src, Path(__file__), Path(__file__).with_name('test_darmanis_contract_20260910.py'),
             scoring_source, MARKERS/'mapping_L1_v3.csv'] + [MARKERS/f'{s}.csv' for s in SETS]
    versions = {}
    for name in ['dgscrna','scanpy','hdbscan','anndata','numpy','scipy','scikit-learn','umap-learn']:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = 'not_available'
    manifest.update({'donor': args.donor, 'slurm_job': os.environ['SLURM_JOB_ID'],
                'dependency_versions':versions, 'scoring_source':str(scoring_source),
                'pilot_selection': '输入细胞数最多的BT_S4；随后固定其余三个donor，不按结果换样本',
                'stage': 'marker_only', 'primary_mcs':15, 'sensitivity_mcs':50,
                'primary_cutoff':'mean', 'cutoffs':['mean','0.5','none'], 'seed':42,
                'pca':30, 'neighbors':15, 'hvg':False, 'top_gene_removal':False,
                'doubletfinder':'not_run_Smart_seq2_no_10x_rate_assumption',
                'deg_top_n':100, 'min_logfc':1., 'max_pval_adj':.05,
                'hashes':{str(p):digest(p) for p in files},
                'mapping_note':'vascular lymphangioblast按已有共同词表归入Vascular；发布标签为concordance reference，不是独立实测gold'})
    (out/'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    sets, mapping = load_panels()
    a = ad.read_h5ad(src)
    a = a[a.obs['donor_id'].astype(str) == args.donor].copy()
    if not a.n_obs:
        raise ValueError('未找到固定donor')
    if not a.obs_names.is_unique:
        raise ValueError('重复cell ID')
    values = a.X.data if sp.issparse(a.X) else np.asarray(a.X).ravel()
    if not (np.isfinite(values).all() and (values >= 0).all() and
            np.isclose(values, np.rint(values), atol=1e-6, rtol=0).all()):
        raise ValueError('X不是非负有限counts-like矩阵；须先核查源数据，禁止继续归一化')
    a.var['source_id'] = a.var_names.astype(str)
    symbols = a.var['feature_name'].astype(str)
    if symbols.duplicated().any() or symbols.isin(['', 'nan']).any():
        raise ValueError('重复/缺失symbol须先审查，不能静默改名导致marker缺失')
    a.var_names = symbols.values
    a.var['mt'] = a.var_names.str.upper().str.startswith('MT-')
    has_mt = bool(a.var['mt'].any())
    sc.pp.calculate_qc_metrics(a, qc_vars=['mt'] if has_mt else [], percent_top=None, log1p=False, inplace=True)
    if not has_mt:
        a.obs['pct_counts_mt'] = np.nan  # 不把未测量伪装成0%。
    manifest['mt_qc'] = 'computed' if has_mt else 'unavailable_in_released_expression_matrix'
    manifest['qc_policy_note'] = 'author_qc使用作者已通过QC的细胞并加nFeature/gene检出过滤，不声称重新完成MT及doublet QC'
    qc = a.obs[['donor_id','Location','Selection','cell_type','total_counts','n_genes_by_counts','pct_counts_mt']].copy()
    qc['retained'] = qc_mask(qc, has_mt, args.qc_policy)
    qc.to_csv(out/'qc_cells.csv')
    manifest['input_shape'] = list(a.shape)
    a = a[qc['retained'].values].copy()
    if a.n_obs <= 31:
        raise ValueError('QC后样本过小，不改变预设阈值抢救')
    sc.pp.filter_genes(a, min_cells=3)
    manifest['qc_shape'] = list(a.shape)
    manifest['counts_provenance'] = 'CELLxGENE下载X，全值counts-like验证通过；不是原始全部测序细胞'
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.layers['lognorm'] = a.X.copy()
    space = a.copy()
    sc.pp.scale(space, max_value=10)
    sc.tl.pca(space, n_comps=30, svd_solver='arpack', random_state=42, mask_var=None)
    sc.pp.neighbors(space, n_neighbors=15, use_rep='X_pca', random_state=42)
    sc.tl.umap(space, random_state=42)
    np.savez_compressed(out/'embedding.npz', pca=space.obsm['X_pca'], umap=space.obsm['X_umap'], cells=np.asarray(a.obs_names,dtype=str))
    gold = a.obs['cell_type'].map(GOLD)
    if gold.isna().any():
        raise ValueError('未定义的作者标签')
    rows, perclass, coverage = [], [], []
    for name, panels in sets.items():
        for panel, genes in panels.items():
            coverage.append({'marker_set':name,'panel':panel,'n_panel_genes':len(genes),
                             'n_detected_unique':len(set(genes)&set(a.var_names)),
                             'target':map_call(name,panel,mapping)})
    pd.DataFrame(coverage).to_csv(out/'panel_coverage.csv', index=False)
    for mcs in (15,50):
        label = hdbscan.HDBSCAN(min_cluster_size=mcs, min_samples=mcs).fit_predict(space.obsm['X_umap'])
        cl = np.asarray(['Noise' if i == -1 else str(i) for i in label])
        a.obs['cluster'] = pd.Categorical(cl)
        groups = sorted(set(cl)-{'Noise'})
        predictions = pd.DataFrame({'cell':a.obs_names,'gold':gold.values,'cluster':cl})
        deg_available = len(groups) >= 2
        if deg_available:
            sc.tl.rank_genes_groups(a, 'cluster', groups=groups, reference='rest', method='wilcoxon',
                                   n_genes=min(100,a.n_vars), layer='lognorm', use_raw=False)
            sc.get.rank_genes_groups_df(a, group=None).to_csv(out/f'deg_mcs{mcs}.csv', index=False)
        for name, panels in sets.items():
            if deg_available:
                score_cell_types(a, panels, 'cluster', marker_set_name=name,
                                 cutoffs=['mean','0.5','none'], min_logfc=1., max_pval_adj=.05)
            reachable = {map_call(name,p,mapping) for p in panels} & set(VOCAB)
            for cutoff in ['mean','0.5','none']:
                key = f'{name}_cluster_{cutoff}'
                raw = a.obs[key].astype(str).values if deg_available else np.repeat('Unknown',a.n_obs)
                pred = np.asarray([map_call(name,p,mapping) for p in raw])
                predictions[f'{name}_{cutoff}_panel'] = raw
                predictions[f'{name}_{cutoff}_label'] = pred
                metrics = evaluate(gold.values,pred)
                metrics.update(marker_set=name,mcs=mcs,cutoff=cutoff,n_clusters=len(groups),
                               status='scored' if deg_available else 'insufficient_clusters',
                               vocabulary_fixed7_upper_bound=len(reachable)/7)
                rows.append(metrics)
                precision, recall, f1, support = precision_recall_fscore_support(gold.values,pred,labels=list(VOCAB),zero_division=0)
                for i,c in enumerate(VOCAB):
                    perclass.append({'marker_set':name,'mcs':mcs,'cutoff':cutoff,'class':c,
                                     'precision':precision[i],'recall':recall[i],'f1':f1[i],'support':int(support[i])})
        predictions.to_csv(out/f'predictions_mcs{mcs}.csv', index=False)
    pd.DataFrame(rows).to_csv(out/'metrics.csv',index=False)
    pd.DataFrame(perclass).to_csv(out/'perclass.csv',index=False)
    print(pd.DataFrame(rows).query("mcs == 15 and cutoff == 'mean'").to_string(index=False),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--donor', choices=DONORS, default='BT_S4')
    parser.add_argument('--run-id', default='pilot_v1')
    parser.add_argument('--qc-policy', choices=['strict_mt','author_qc'], default='strict_mt')
    run(parser.parse_args())
