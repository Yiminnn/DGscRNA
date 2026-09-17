#!/usr/bin/env python3
"""CM2_glioma_other bounded cluster x cutoff grid on the GBM-CARE cohort.

Protocol: results/g274_cgo_grid_v1/protocol.json. Features, marker set and scoring are
frozen; only the clusterer and the abstention cutoff move. Two deliberate reuses keep the
numbers on the existing ruler instead of a new one:

  * scoring  -- handoff/g274_table4/v5_final_annotations.py (label_scope/metrics), which
                itself AST-extracts the legacy Lfine n>=20 helpers from handoff/g274/grid_sets.py
  * DL stage -- handoff/gbm/darmanis_region_final.py refine_final, the corrected semantics
                (Noise never trains and stays terminal Unknown, non-noise abstentions are
                normalised to Undecided, failures are flagged instead of silently scored)

The feature path is copied from handoff/g274/annotate_v3.py, which produced the frozen
baseline; that script cannot be imported (it reads sys.argv at module level).
"""
from runtime_guard import bootstrap

if __name__ == '__main__':
    bootstrap()

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT = ROOT / 'results/g274_cgo_grid_v1'
MARKERS = ROOT / 'handoff/markers_v3'
MTX = ROOT / 'data_bench/GSE274546/mtx'
ROSTER = ROOT / 'results/g274_v5/final_annotations/cohort_verified_final_samples.csv'
MARKER = 'CM2_glioma_other'
SEED, NOISE = 42, 'Noise'
CUTOFFS = ('none', 'mean', '0.5')
CONFIGS = ('umap_hdb15', 'pca_hdb15', 'umap_hdb5',
           'leiden_0p5', 'leiden_1p0', 'leiden_2p0', 'pca_kmeans7')
sys.path.insert(0, str(ROOT / 'handoff/gbm'))


def roster():
    return pd.read_csv(ROSTER)['sample'].tolist()


def digest(path):
    with open(path, 'rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


# ---------------------------------------------------------------- scoring ruler

def load_helpers():
    import v5_final_annotations as v5
    helpers = v5.legacy_helpers()
    mapping = pd.read_csv(MARKERS / 'mapping_L1_v3.csv')
    helpers['panel_labels'] = {MARKER: set(mapping.loc[mapping.marker_set.eq(MARKER), 'gold_class'])}
    return helpers


def reference_frame(sample):
    """Per-cell reference columns in input order; Lfine follows annotate_v3.py."""
    obs = pd.read_csv(MTX / sample / 'obs.csv')
    fine = np.where(obs.L1.astype(str) == 'Malignant',
                    'Malignant_' + obs.MalState.astype(str), obs.L3.astype(str))
    fine = pd.Series(fine, index=obs.index).replace('nan', np.nan)
    fine = fine.fillna(obs.L2.astype(str).replace('nan', np.nan)).fillna(obs.L1.astype(str))
    obs = obs.copy()
    obs['Lfine'] = fine.values
    return obs


def label_scope(pc, helpers):
    import v5_final_annotations as v5
    return v5.label_scope(pc, helpers)


def score_prediction(pred, pc, scope, helpers):
    """Legacy Lfine set-valued macro-F1 and its companion columns, unchanged."""
    import v5_final_annotations as v5
    pred = np.asarray([str(p) for p in pred], dtype=object)
    assert len(pred) == len(pc), 'prediction lost cells'
    return v5.metrics(pred, pc, scope, helpers, MARKER)


# ---------------------------------------------------------------- pipeline

def build_adata(sample):
    """annotate_v3.py feature path: HVG2000 -> scale -> PCA30 -> kNN15 -> UMAP2."""
    import anndata as ad
    import scanpy as sc
    import scipy.io as sio

    src = MTX / sample
    X = sio.mmread(str(src / 'matrix.mtx')).T.tocsr()
    genes = [line.strip() for line in open(src / 'genes.tsv')]
    obs = pd.read_csv(src / 'obs.csv')
    a = ad.AnnData(X.astype(np.float32))
    a.var_names = genes
    a.obs = obs.set_index('CellID')
    a.var_names_make_unique()
    sc.pp.filter_genes(a, min_cells=3)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.layers['lognorm'] = a.X.copy()
    sc.pp.highly_variable_genes(a, n_top_genes=2000)
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=30, random_state=SEED)
    # annotate_v3 re-scales a copy of the same lognorm matrix for the DL features; reuse the
    # matrix already in .X instead of copying the whole object -- identical numbers, less RAM.
    a.obsm['X_scaled'] = np.asarray(a.X, dtype=np.float32)
    a.X = a.layers['lognorm'].copy()
    sc.pp.neighbors(a, n_neighbors=15, random_state=SEED)
    sc.tl.umap(a, random_state=SEED)
    return a


def cluster_labels(a, config):
    import hdbscan
    import scanpy as sc
    from sklearn.cluster import KMeans

    if config.startswith('leiden_'):
        resolution = float(config.split('_')[1].replace('p', '.'))
        sc.tl.leiden(a, resolution=resolution, flavor='igraph', n_iterations=2,
                     directed=False, random_state=SEED, key_added='_leiden')
        return np.asarray([str(v) for v in a.obs['_leiden']], dtype=object)
    if config == 'pca_kmeans7':
        raw = KMeans(n_clusters=7, n_init=10, random_state=SEED).fit_predict(a.obsm['X_pca'])
    else:
        space = a.obsm['X_umap'] if config.startswith('umap_') else a.obsm['X_pca']
        size = 5 if config == 'umap_hdb5' else 15
        raw = hdbscan.HDBSCAN(min_cluster_size=size, min_samples=size).fit_predict(space)
    return np.array([NOISE if c == -1 else str(c) for c in raw], dtype=object)


def load_panels(a):
    """Panel genes restricted to the measured vocabulary, as annotate_v3.py does."""
    table = pd.read_csv(MARKERS / f'{MARKER}.csv', index_col=0)
    present = set(a.var_names)
    return {c: [g for g in table[c].dropna().tolist() if g in present] for c in table.columns}


def score_panels(deg, panels, order):
    """R density_score: sum of logFC over panel genes in the cluster DE list / panel length."""
    names = list(panels)
    S = np.zeros((len(names), len(order)), dtype=np.float32)
    for j, cluster in enumerate(order):
        d = deg[cluster]
        for i, panel in enumerate(names):
            genes = panels[panel]
            S[i, j] = (sum(d[g] for g in genes if g in d) / len(genes)) if genes else 0.0
    return S, names


def calls_at_cutoff(S, names, cutoff):
    """R source.R:410-418. 'mean' is a threshold over the per-cluster maxima of this sample."""
    threshold = float(S.max(0).mean())
    calls = []
    for j in range(S.shape[1]):
        column = S[:, j]
        top = float(column.max())
        winners = np.flatnonzero(column == top)
        if top <= 0 or len(winners) > 1:
            calls.append('Undecided')
            continue
        call = names[int(winners[0])]
        if cutoff == 'mean' and top < threshold:
            call = 'Undecided'
        elif cutoff == '0.5' and top < 0.5:
            call = 'Undecided'
        calls.append(call)
    return calls


def seed_labels(cl, call_by_cluster, mapping):
    """Cluster calls -> per-cell gold labels. UNMAPPABLE is an abstention, not a class."""
    projected = {c: ('Undecided' if v == 'Undecided' else mapping.get(v, 'UNMAPPABLE'))
                 for c, v in call_by_cluster.items()}
    seed = np.array([projected.get(c, 'Unknown') for c in cl], dtype=object)
    return np.where(seed == 'UNMAPPABLE', 'Undecided', seed).astype(object)


def final_prediction(a, cl, seed, train_fn=None, predict_fn=None):
    """Corrected DL stage: Noise terminal, abstentions pooled, failures flagged."""
    import anndata as ad
    import scipy.sparse as sp
    from darmanis_region_final import refine_final

    # The classifier only ever reads obsm['X_scaled']; an empty sparse .X of the right width
    # keeps refine_final's feature-width assertion honest without a second dense copy.
    slim = ad.AnnData(sp.csr_matrix((a.n_obs, a.n_vars), dtype=np.float32))
    slim.obs_names = a.obs_names
    slim.var_names = a.var_names
    slim.obs['cluster'] = cl
    slim.obs['seed_panel'] = seed
    slim.obsm['X_scaled'] = a.obsm['X_scaled']
    return refine_final(slim, 'seed_panel', train_fn, predict_fn)


# ---------------------------------------------------------------- driver

def run_sample(sample, outdir):
    import torch
    torch.set_num_threads(1)

    started = time.time()
    helpers = load_helpers()
    mapping = pd.read_csv(MARKERS / 'mapping_L1_v3.csv')
    mapping = mapping[mapping.marker_set.eq(MARKER)].set_index('panel')['gold_class'].to_dict()
    pc = reference_frame(sample)
    scope = label_scope(pc, helpers)
    a = build_adata(sample)
    assert a.n_obs == len(pc) and list(a.obs_names) == list(pc.CellID), 'cell roster drift'
    panels = load_panels(a)

    rows, cells = [], {'CellID': list(pc.CellID)}
    for config in CONFIGS:
        elapsed = time.time()
        try:
            cl = cluster_labels(a, config)
            order = [c for c in pd.Series(cl).unique().tolist() if c != NOISE]
            cells[f'cluster_{config}'] = cl
            if len(order) < 2:
                status = 'all_noise' if not order else 'lt2_clusters'
                for cutoff in CUTOFFS:
                    pred = np.array(['Unknown'] * a.n_obs, dtype=object)
                    cells[f'pred_{config}__{cutoff}'] = pred
                    rows.append(dict(sample=sample, config=config, cutoff=cutoff,
                                     cluster_status=status, n_clusters=len(order),
                                     n_noise=int((cl == NOISE).sum()), dl_status='not_reached',
                                     final_valid=False, n_pool=0, n_training=0,
                                     n_training_classes=0, n_dl_assigned=0,
                                     seconds=round(time.time() - elapsed, 1),
                                     **score_prediction(pred, pc, scope, helpers)))
                continue
            import scanpy as sc
            a.obs['_cl'] = pd.Categorical(cl)
            sc.tl.rank_genes_groups(a, '_cl', groups=order, method='wilcoxon', n_genes=100,
                                    layer='lognorm', use_raw=False)
            deg = {c: dict(zip(a.uns['rank_genes_groups']['names'][c],
                               a.uns['rank_genes_groups']['logfoldchanges'][c])) for c in order}
            S, names = score_panels(deg, panels, order)
            np.save(outdir / f'scorematrix_{config}.npy', S)
            for cutoff in CUTOFFS:
                seed = seed_labels(cl, dict(zip(order, calls_at_cutoff(S, names, cutoff))), mapping)
                final, _, info = final_prediction(a, cl, seed)
                pred = np.asarray(final.astype(str).to_numpy(), dtype=object)
                cells[f'pred_{config}__{cutoff}'] = pred
                rows.append(dict(sample=sample, config=config, cutoff=cutoff,
                                 cluster_status='ok', n_clusters=len(order),
                                 n_noise=int((cl == NOISE).sum()),
                                 dl_status=info['dl_status'], final_valid=bool(info['final_valid']),
                                 n_pool=info['n_pool'], n_training=info['n_training'],
                                 n_training_classes=info['n_training_classes'],
                                 n_dl_assigned=info['n_dl_assigned'],
                                 dl_error=info.get('error', ''),
                                 seconds=round(time.time() - elapsed, 1),
                                 **score_prediction(pred, pc, scope, helpers)))
                print(f'  {sample} {config} {cutoff} pool={info["n_pool"]} '
                      f'{info["dl_status"]} F1={rows[-1]["lfine_macroF1"]:.4f} '
                      f'{time.time() - started:.0f}s', flush=True)
        except Exception as exc:  # a broken config must not take the sample down
            note = f'error:{type(exc).__name__}: {str(exc)[:120]}'
            print(f'  {sample} {config} {note}\n{traceback.format_exc()}', flush=True)
            for cutoff in CUTOFFS:
                pred = np.array(['Unknown'] * a.n_obs, dtype=object)
                cells.setdefault(f'pred_{config}__{cutoff}', pred)
                rows.append(dict(sample=sample, config=config, cutoff=cutoff,
                                 cluster_status=note, n_clusters=np.nan, n_noise=np.nan,
                                 dl_status='not_reached', final_valid=False,
                                 seconds=round(time.time() - elapsed, 1),
                                 **score_prediction(pred, pc, scope, helpers)))

    table = pd.DataFrame(rows)
    assert len(table) == len(CONFIGS) * len(CUTOFFS), 'grid incomplete'
    table.to_csv(outdir / 'grid.csv', index=False)
    pd.DataFrame(cells).to_parquet(outdir / 'percell.parquet', index=False)
    np.save(outdir / 'umap.npy', a.obsm['X_umap'])
    manifest = dict(sample=sample, n_cells=int(a.n_obs), n_genes=int(a.n_vars),
                    lfine_n_classes=len(scope[1]), marker_set=MARKER, seed=SEED,
                    slurm_job_id=os.environ.get('SLURM_JOB_ID'), host=os.uname().nodename,
                    seconds=round(time.time() - started, 1),
                    sources={str(p): digest(p) for p in
                             [MTX / sample / 'matrix.mtx', MTX / sample / 'obs.csv',
                              MTX / sample / 'genes.tsv', MARKERS / f'{MARKER}.csv',
                              MARKERS / 'mapping_L1_v3.csv', Path(__file__).resolve(),
                              ROOT / 'handoff/g274/grid_sets.py',
                              ROOT / 'handoff/g274_table4/v5_final_annotations.py',
                              ROOT / 'handoff/gbm/darmanis_region_final.py']})
    (outdir / 'manifest.json').write_text(json.dumps(manifest, indent=2, default=str) + '\n')
    print(f'{sample}: DONE {len(table)} arms in {time.time() - started:.0f}s', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample')
    parser.add_argument('--index', type=int, help='1-based row of the frozen roster')
    args = parser.parse_args()
    assert os.environ.get('SLURM_JOB_ID'), 'scientific compute runs on SLURM only'
    sample = args.sample or roster()[args.index - 1]
    assert sample in roster(), f'{sample} is not on the frozen 111-sample roster'
    destination = OUT / 'per_sample' / sample
    staging = OUT / 'per_sample' / f'.tmp_{sample}_{os.environ["SLURM_JOB_ID"]}'
    shutil.rmtree(staging, ignore_errors=True)
    staging.mkdir(parents=True)
    run_sample(sample, staging)
    shutil.rmtree(destination, ignore_errors=True)
    staging.rename(destination)
    print(f'PLACED {destination}', flush=True)


if __name__ == '__main__':
    main()
