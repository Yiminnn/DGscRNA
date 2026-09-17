"""Plot every saved GBM partition and terminal annotation without refitting."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
BASE = ROOT / 'results/hvg_ptc_20260916_v1'
DEST = BASE / 'gbm_cluster_figures'
FEATURES = ['all', 'hvg500', 'hvg1000', 'hvg2000', 'hvg3000', 'hvg5000']
REDUCERS = ['PCA', 'FA', 'ICA', 'Isomap', 'TSNE', 'UMAP', 'none']
METHODS = ['KMeans', 'GMM', 'HDBSCAN']
CL_COLORS = ['#0072B2', '#E69F00', '#009E73', '#CC79A7', '#56B4E9', '#D55E00', '#000000', '#999933']
L1_COLORS = {
    'Malignant': '#332288', 'TAM': '#117733', 'Lymphocyte': '#44AA99',
    'Oligodendrocyte': '#88CCEE', 'Astrocyte': '#DDCC77', 'OPC': '#CC6677',
    'Excitatory neuron': '#AA4499', 'Inhibitory neuron': '#882255',
    'Endothel': '#661100', 'Pericyte': '#6699CC', 'Other': '#888888',
    'Unknown': '#C5C5C5', 'Undecided': '#C5C5C5', 'Noise': '#C5C5C5'}


def sha(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def save(p, obj):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + '.part')
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    tmp.replace(p)


def arm_label(a):
    if a['clusterer'] == 'HDBSCAN':
        return f"HDBSCAN {a['min_cluster_size']}/{a['min_samples']}"
    return f"{a['clusterer']} K{a['k']}" + (f" {a['covariance']}" if a['clusterer'] == 'GMM' else '')


def geometry_label(g):
    if g['dr'] == 'none':
        return 'No DR [PCA2 display only]'
    prefix = 'PCA30 > ' if g['input_space'] == 'pca30' else ''
    return prefix + g['dr'] + str(g['dim']) + (' [axes 1-2]' if g['dim'] > 2 else '')


def main(sample, limit=None):
    assert os.environ.get('SLURM_JOB_ID'), 'All plotting and array reads require SLURM'
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.lines import Line2D
    import re
    import subprocess

    started = time.time()
    out = DEST / 'samples' / sample
    if limit:
        out = DEST / 'preview' / (sample + '_' + sha(__file__)[:8])
    out.mkdir(parents=True, exist_ok=True)
    sources = {str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__), BASE / 'protocol/geometries.json',
        BASE / 'evaluation' / sample / 'metrics.csv', BASE / 'prepared' / sample / 'cells.tsv']}
    if (out / 'COMPLETE').exists():
        m = json.loads((out / 'manifest.json').read_text())
        assert (out / 'COMPLETE').read_text().strip() == sha(out / 'manifest.json')
        assert m['sources'] == sources
        assert all(sha(out / name) == value for name, value in m['outputs'].items())
        print(sample, 'verified existing figures', flush=True)
        return

    gs = json.loads((BASE / 'protocol/geometries.json').read_text())
    by_geo = {g['geometry_id']: g for g in gs}
    rows = [(g, a) for g in gs for a in g['arms']]
    rows.sort(key=lambda ga: ('E1_primary_table' not in ga[1]['families'],
        FEATURES.index(ga[0]['feature']) if ga[0]['feature'] in FEATURES else 10,
        ga[0]['feature'], REDUCERS.index(ga[0]['dr']) if ga[0]['dr'] in REDUCERS else 10,
        ga[0]['input_space'], ga[0]['dim'] or 0, ga[0]['seed'],
        METHODS.index(ga[1]['clusterer']), ga[1]['arm_id'], ga[0]['geometry_id']))
    assert len(rows) == 330 and sum('E1_primary_table' in a['families'] for g, a in rows) == 126
    cells_path = BASE / 'prepared' / sample / 'cells.tsv'
    cells = cells_path.read_text().splitlines()
    raw_obs = ROOT / 'data_bench/GSE274546/mtx' / sample / 'obs.csv'
    frozen_input = json.loads((ROOT / 'results/g274_cohort/inputs_v1' / sample / 'manifest.json').read_text())
    assert sha(raw_obs) == frozen_input['sha256'][str(raw_obs)]
    obs = pd.read_csv(raw_obs).set_index('CellID')
    assert obs.index.is_unique and set(cells) <= set(obs.index)
    ref = obs.loc[cells, 'L1'].to_numpy(dtype=str)
    assert set(ref) <= set(L1_COLORS)
    metric = pd.read_csv(BASE / 'evaluation' / sample / 'metrics.csv').set_index('condition')
    assert len(metric) == 330 and metric.index.is_unique
    evaluation_manifest = json.loads((BASE / 'evaluation' / sample / 'manifest.json').read_text())
    assert evaluation_manifest['outputs']['metrics.csv'] == sha(BASE / 'evaluation' / sample / 'metrics.csv')
    assert all(metric.n_cells == len(cells))
    pca_g = next(g for g in gs if g['feature'] == 'all' and g['dr'] == 'PCA' and g['dim'] == 2 and g['seed'] == 42)
    embeddings, data, filehashes = {}, {}, {}

    def remember(p, expected=None):
        h = sha(p)
        assert expected is None or expected == h, str(p)
        filehashes[str(p.relative_to(ROOT))] = h
        return h

    def xy(g):
        real = pca_g if g['dr'] == 'none' else g
        gid = real['geometry_id']
        if gid not in embeddings:
            folder = BASE / 'fits' / sample / gid
            em = folder / 'embedding_manifest.json'
            path = folder / 'embedding.npy'
            if not path.exists():
                return None
            emeta = json.loads(em.read_text())
            remember(em)
            remember(path, emeta['sha256'])
            value = np.load(path, allow_pickle=False)
            assert value.ndim == 2 and value.shape[0] == len(cells) and value.shape[1] >= 2
            assert np.isfinite(value).all()
            embeddings[gid] = value[:, :2].copy()
        return embeddings[gid]

    for g, a in rows:
        cid = g['geometry_id'] + '/' + a['arm_id']
        folder = BASE / 'fits' / sample / cid
        row = metric.loc[cid]
        record = dict(g=g, a=a, cid=cid, metric=row, xy=xy(g), cl=None, final=None,
                      cluster_source=None, final_source=None)
        if (folder / 'clusters.npy').exists():
            am_path = folder / 'manifest.json'
            am = json.loads(am_path.read_text())
            remember(am_path)
            remember(folder / 'clusters.npy', am.get('outputs', {}).get('clusters.npy'))
            cl = np.load(folder / 'clusters.npy', allow_pickle=False)
            assert cl.shape == (len(cells),) and np.issubdtype(cl.dtype, np.integer)
            assert record['xy'] is not None
            record.update(cl=cl, cluster_source=str((folder / 'clusters.npy').relative_to(ROOT)))
        if row['status'] == 'completed':
            assert (folder / 'COMPLETE').read_text().strip() == remember(folder / 'manifest.json')
            assert evaluation_manifest['result_manifest_hashes'][cid] == sha(folder / 'manifest.json')
            remember(folder / 'predictions.npz', am['outputs']['predictions.npz'])
            with np.load(folder / 'predictions.npz', allow_pickle=False) as pred:
                np.testing.assert_array_equal(pred['cluster'], record['cl'])
                final = pred['final'].copy()
            assert final.shape == (len(cells),) and set(final) <= set(L1_COLORS)
            assert (final[record['cl'] == -1] == 'Unknown').all()
            record.update(final=final, final_source=str((folder / 'predictions.npz').relative_to(ROOT)))
        else:
            assert str(row['status']).startswith(('structural', 'numerical')), row['status']
        data[cid] = record

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
        'axes.titlesize': 8, 'axes.labelsize': 7, 'xtick.labelsize': 6,
        'ytick.labelsize': 6, 'pdf.fonttype': 42, 'savefig.facecolor': 'white'})
    reference_legend = [Line2D([], [], marker='o', color='none', markerfacecolor=color,
        markeredgewidth=0, markersize=5, label=label) for label, color in L1_COLORS.items()
        if label not in ['Undecided', 'Noise']]

    def panel(ax, rec, endpoint, title=None, labels=True):
        values = rec['cl'] if endpoint == 'clusters' else (ref if endpoint == 'reference' else rec['final'])
        z = rec['xy']
        ax.set_title(title or '', fontsize=7.5, pad=4)
        if values is None or z is None:
            reason = str(rec['metric'].get('structural_reason', ''))
            if reason == 'nan': reason = 'Verified numerical fitting failure'
            ax.text(.5, .53, 'Unavailable', ha='center', va='center', transform=ax.transAxes, fontsize=10)
            import textwrap
            ax.text(.5, .34, '\n'.join(textwrap.wrap(reason, 38)), ha='center', va='center',
                    transform=ax.transAxes, fontsize=6.5, color='#555555')
            ax.set_facecolor('#F2F2F2')
            ax.set_xticks([]); ax.set_yticks([])
            return
        colors = [('#C5C5C5' if int(v) == -1 else CL_COLORS[int(v) % len(CL_COLORS)])
                  for v in values] if endpoint == 'clusters' else [L1_COLORS[v] for v in values]
        ax.scatter(z[:, 0], z[:, 1], c=colors, s=max(.25, min(2.2, 9000 / len(cells))),
                   alpha=.72, linewidths=0, rasterized=True)
        ax.set_aspect('equal', adjustable='datalim')
        # Numeric labels distinguish clusters even when the eight-color palette repeats.
        if endpoint == 'clusters' and labels:
            for label in np.unique(values):
                if label == -1: continue
                center = np.median(z[values == label], axis=0)
                ax.text(*center, str(label), fontsize=6, ha='center', va='center', color='black',
                    path_effects=[pe.withStroke(linewidth=1.7, foreground='white')])
        ax.set_xticks([]); ax.set_yticks([])
        ax.spines[['top', 'right']].set_visible(False)
        ax.set_xlabel('PCA1 (display only)' if rec['g']['dr'] == 'none' else rec['g']['dr'] + '1', labelpad=1)
        ax.set_ylabel('PCA2 (display only)' if rec['g']['dr'] == 'none' else rec['g']['dr'] + '2', labelpad=1)

    ordered = [data[g['geometry_id'] + '/' + a['arm_id']] for g, a in rows]
    if limit: ordered = ordered[:limit]
    index_rows = []
    for i, rec in enumerate(ordered):
        r = rec['metric']
        index_rows.append(dict(sample=sample, condition=rec['cid'],
            geometry_id=rec['g']['geometry_id'], arm_id=rec['a']['arm_id'],
            feature=rec['g']['feature'], representation=geometry_label(rec['g']),
            clusterer=rec['a']['clusterer'], cluster_parameters=arm_label(rec['a']),
            seed=rec['g']['seed'], scoring_features=rec['a']['scoring_features'],
            neighbors=rec['g']['neighbors'], min_dist=rec['g']['min_dist'],
            k=rec['a']['k'], covariance=rec['a']['covariance'],
            min_cluster_size=rec['a']['min_cluster_size'], min_samples=rec['a']['min_samples'],
            families=';'.join(rec['a']['families']), status=r['status'],
            n_cells=len(cells), cluster_plotted=rec['cl'] is not None,
            terminal_plotted=rec['final'] is not None, clusters_pdf='clusters_all_conditions.pdf',
            terminal_pdf='terminal_all_conditions.pdf', pdf_page=i // 12 + 1, panel=i % 12 + 1,
            cluster_source=rec['cluster_source'], terminal_source=rec['final_source'],
            displayed_coordinates='saved all-gene PCA2, display only' if rec['g']['dr'] == 'none'
                else 'first two saved fitted coordinates',
            partition_ari=r.get('partition_ari'), terminal_strict_L1_macroF1_present=r.get('terminal_strict_L1_macroF1_present')))

    for endpoint in ['clusters', 'terminal']:
        path = out / (endpoint + '_all_conditions.pdf')
        with PdfPages(path, metadata={'Title': f'{sample}: all GBM {endpoint}',
              'Subject': 'Saved fits; every condition and every cell retained; no refitting'}) as pdf:
            for page, start in enumerate(range(0, len(ordered), 12), 1):
                selected = ordered[start:start + 12]
                fig, axs = plt.subplots(4, 3, figsize=(11.7, 11.2))
                fig.subplots_adjust(left=.045, right=.99, bottom=.08, top=.91, hspace=.6, wspace=.15)
                for ax, rec in zip(axs.flat, selected):
                    title = (f"{rec['g']['feature']} | {geometry_label(rec['g'])}\n"
                        f"{arm_label(rec['a'])} | seed {rec['g']['seed']} | scoring {rec['a']['scoring_features']}")
                    if rec['g']['dr'] in ['UMAP', 'SCANPY_UMAP', 'Isomap']:
                        title += f"\nneighbors={rec['g']['neighbors']}"
                        if 'UMAP' in rec['g']['dr']:
                            title += f" | min_dist={rec['g']['min_dist']}"
                    panel(ax, rec, endpoint, title)
                for ax in list(axs.flat)[len(selected):]: ax.axis('off')
                title = 'Cluster IDs (partition diagnostic)' if endpoint == 'clusters' else 'Terminal DG-scRNA DL/refinement annotation'
                fig.suptitle(f'{sample} | {title}\nAll {len(cells):,} cells; conditions {start + 1}-{start + len(selected)} / {len(ordered)}; page {page}', y=.975, fontsize=12)
                footer = ('Gray = HDBSCAN noise; numeric cluster IDs are local to each fit. Repeated colors are distinguished by IDs.\n'
                    if endpoint == 'clusters' else 'Gray = final Unknown/abstention. Missing terminal outputs are explicitly unavailable, never replaced by marker-only calls.\n')
                footer += 'Native 2D fits: saved coordinates. Fits above 2D: axes 1-2 only. No DR: fixed all-gene PCA2 display; fitting used all selected dimensions.'
                fig.text(.5, .015, footer, ha='center', va='bottom', fontsize=6.5)
                if endpoint == 'terminal':
                    fig.legend(handles=reference_legend, loc='lower center', bbox_to_anchor=(.5, .041), ncol=7, frameon=False, fontsize=6.5)
                pdf.savefig(fig, dpi=150)
                # Inline notebook supplements: every non-primary condition for the fixed example.
                if sample == 'TKU3186' and endpoint == 'clusters' and start + len(selected) > 126:
                    fig.savefig(out / f'secondary_clusters_page_{page:02d}.png', dpi=135)
                if page == 1:
                    fig.savefig(out / f'{endpoint}_preview.png', dpi=110)
                plt.close(fig)
        info = subprocess.check_output(['pdfinfo', str(path)], text=True)
        assert int(re.search(r'^Pages:\s+(\d+)', info, re.M).group(1)) == (len(ordered) + 11) // 12
        pages = subprocess.check_output(['pdftotext', str(path), '-'], text=True).split('\f')
        assert all(sample in page for page in pages if page.strip())
        print(sample, endpoint, 'atlas complete', flush=True)

    primary_images = []
    if sample == 'TKU3186' and not limit:
        for feature in FEATURES:
            for endpoint in ['clusters', 'terminal']:
                fig, axs = plt.subplots(7, 4, figsize=(13, 17.8))
                fig.subplots_adjust(left=.055, right=.992, top=.93, bottom=.065, hspace=.30, wspace=.10)
                for ri, reducer in enumerate(REDUCERS):
                    selected = [r for r in ordered if r['g']['feature'] == feature and r['g']['dr'] == reducer
                                and 'E1_primary_table' in r['a']['families']]
                    assert len(selected) == 3
                    selected.sort(key=lambda r: METHODS.index(r['a']['clusterer']))
                    panel(axs[ri, 0], selected[0], 'reference', 'Saved L1 reference' if ri == 0 else '')
                    for ci, rec in enumerate(selected, 1):
                        score_key = 'partition_ari' if endpoint == 'clusters' else 'terminal_strict_L1_macroF1_present'
                        score = rec['metric'].get(score_key)
                        scoretext = 'unavailable' if pd.isna(score) else f'{score:.3f}'
                        short = 'ARI' if endpoint == 'clusters' else 'L1 F1'
                        title = f"{METHODS[ci-1]} | {short} {scoretext}"
                        panel(axs[ri, ci], rec, endpoint, title)
                    axs[ri, 0].text(-.18, .5, reducer if reducer != 'none' else 'No DR', transform=axs[ri, 0].transAxes,
                        rotation=90, ha='center', va='center', fontsize=11, fontweight='bold')
                title = 'cluster partitions' if endpoint == 'clusters' else 'terminal DL/refinement calls'
                fig.suptitle(f'TKU3186 | {feature} | all 7 representations x 3 clusterers\n{title}; every one of {len(cells):,} cells is shown', fontsize=15, y=.98)
                fig.legend(handles=reference_legend, loc='lower center', bbox_to_anchor=(.5, .028), ncol=7, frameon=False, fontsize=8)
                fig.text(.5, .009, 'Reference labels are evaluation annotations, not independent truth. Every row uses one frozen geometry.\n'
                    'ARI uses the archived fine labels; F1 uses strict L1. No DR: full selected expression for fitting, PCA2 only for display.\n'
                    'Cluster IDs are fit-specific; gray indicates noise/Unknown.',
                    ha='center', fontsize=8)
                stem = f'primary_{feature}_{endpoint}'
                fig.savefig(out / (stem + '.pdf'), dpi=200)
                fig.savefig(out / (stem + '.png'), dpi=145)
                plt.close(fig)
                primary_images.append(stem)
    pd.DataFrame(index_rows).to_csv(out / 'condition_figure_index.csv', index=False)
    save(out / 'input_hashes.json', filehashes)
    output_hashes = {p.name: sha(p) for p in out.iterdir() if p.is_file() and p.name not in ['manifest.json', 'COMPLETE']}
    m = dict(status='complete', sample=sample, n_conditions=len(ordered), n_cells=len(cells),
        n_cluster_plots=sum(r['cl'] is not None for r in ordered),
        n_terminal_plots=sum(r['final'] is not None for r in ordered),
        atlas_pages=(len(ordered)+11)//12, primary_images=primary_images,
        sources=sources, outputs=output_hashes, slurm_job=os.environ['SLURM_JOB_ID'],
        plotted_all_cells=True, fit_and_annotation_sources_read_only=True,
        elapsed_seconds=time.time()-started)
    save(out / 'manifest.json', m)
    (out / 'COMPLETE').write_text(sha(out / 'manifest.json') + '\n')
    print(json.dumps({k:v for k,v in m.items() if k not in ['sources', 'outputs']}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample')
    parser.add_argument('--index', type=int)
    parser.add_argument('--limit', type=int)
    args = parser.parse_args()
    samples = (BASE / 'protocol/samples.txt').read_text().split()
    index = args.index if args.index is not None else int(os.environ.get('SLURM_ARRAY_TASK_ID', '0'))
    main(args.sample or samples[index], args.limit)
