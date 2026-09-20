"""Evaluate only frozen predictions, and plot every candidate on fixed coordinates."""
from pathlib import Path
import json
import os
import shutil
import sys
CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE / 'legacy'))
from common import OUT as REFERENCE, L1, require_slurm, write_json, sha, checked, complete, utc


def run(directory):
    require_slurm()
    import numpy as np
    import pandas as pd
    import evaluate
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    directory = Path(directory)
    assert checked(directory, 'fit_manifest.json', 'FIT_COMPLETE')
    cfg = json.loads((directory / 'config.json').read_text())
    sample, budget, space = cfg['sample'], cfg['budget'], cfg['space']
    native_cells = Path(cfg['prep']) / 'cells.csv'
    if not (directory / 'cells.csv').exists():
        shutil.copy2(native_cells, directory / 'cells.csv')
    assert sha(directory / 'cells.csv') == sha(native_cells)
    if checked(directory / 'evaluation') and checked(directory / 'figures'):
        return
    write_json(directory / 'prepare_manifest.json', dict(sample=sample, budget=budget, seed=42,
        provenance='Evaluation wrapper only; original native RNA preparation retained at cfg.prep',
        original_prepare_manifest_sha256=sha(Path(cfg['prep']) / 'prepare_manifest.json')))
    evaluate.ROUTES = [Path(c['dest']).name for c in cfg['conditions']]
    evaluate.run(directory, only_arms=['L00_mean'])
    metrics = pd.read_csv(directory / 'evaluation/metrics.csv', dtype={'cutoff': str})
    assert len(metrics) == 3 * len(cfg['conditions'])
    display = REFERENCE / 'GBM' / sample / 'hvg2000/UMAP2.csv'
    frame = pd.read_csv(display, index_col=0)
    coords = frame.to_numpy()
    truth = pd.read_csv(REFERENCE / 'evaluation_inputs' / sample / 'truth.csv.gz', dtype=str, keep_default_na=False)
    assert list(frame.index) == list(truth.cell_id)
    mapping = pd.read_csv(REFERENCE / 'markers/panel_L1_mapping.csv', dtype=str, keep_default_na=False)
    mapping = mapping[mapping.library == 'CM2_glioma_other']
    lookup = dict(zip(mapping.panel, mapping.L1))
    cm = plt.get_cmap('tab20')
    colors = {label: cm(i) for i, label in enumerate(L1)}
    colors.update(Unknown='#bdbdbd', UNMAPPABLE='#7b614c', AMBIGUOUS_NEURON='#8c8c33', NO_L1_COUNTERPART='#7b614c')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'pdf.fonttype': 42})

    def panel(ax, labels, title, cluster=False, hdb=False):
        labels = np.asarray(labels, dtype=str)
        values = sorted(set(labels), key=int if cluster else None)
        palette = {label: cm(i % 20) for i, label in enumerate(values)} if cluster else colors
        if cluster and hdb:
            palette['0'] = '#bdbdbd'
        for label in values:
            take = labels == label
            ax.scatter(coords[take, 0], coords[take, 1], s=2.8, c=[palette.get(label, '#7b614c')], linewidths=0, rasterized=True, alpha=.8)
            if cluster:
                ax.text(*np.median(coords[take], axis=0), 'noise' if hdb and label == '0' else label,
                        fontsize=6, ha='center', bbox=dict(facecolor='white', alpha=.7, edgecolor='none', pad=.5))
        ax.set_title(title, fontsize=9)
        ax.set_xticks([]);ax.set_yticks([])
        ax.set_xlabel('Fixed original HVG2000 UMAP1', fontsize=7);ax.set_ylabel('UMAP2', fontsize=7)
        for spine in ax.spines.values():
            spine.set_visible(False)

    figdir = directory / 'figures';figdir.mkdir(exist_ok=True)
    records = []
    for condition in cfg['conditions']:
        source = Path(condition['dest']);name = source.name
        tm = json.loads((source / 'terminal/L00_mean/terminal_manifest.json').read_text())
        cl = pd.read_csv(source / 'clusters.csv', dtype=str, keep_default_na=False)
        pp = pd.read_csv(source / 'terminal/L00_mean/predictions.csv.gz', dtype=str, keep_default_na=False)
        assert list(cl.cell_id) == list(pp.cell_id) == list(truth.cell_id)
        pred = [lookup.get(v, 'Unknown' if v in ['Unknown', 'Undecided', 'Noise', ''] else 'UNMAPPABLE') for v in pp.final090]
        row = metrics[(metrics.route == name) & (metrics.stage == 'terminal090')].iloc[0]
        fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.3), layout='constrained')
        panel(axes[0], truth.L1, 'Original author L1')
        panel(axes[1], cl.cluster, f'{space} / {name}: partition', True, condition['method'] == 'HDBSCAN_R')
        stage_title = 'Final DL threshold 0.90' if tm['training_executed'] else 'Terminal result: DL not executed'
        panel(axes[2], pred, f'{stage_title}\n{tm["dl_status"]}\nmacro-F1 {row.macroF1_present:.3f}; coverage {row.coverage:.1%}')
        fig.suptitle(f'{sample} | {budget} | {space} / {name}\nCM2_glioma_other / mean; display is shared, not the fitting input', fontsize=11)
        labels = [v for v in L1 + ['Unknown', 'UNMAPPABLE', 'AMBIGUOUS_NEURON', 'NO_L1_COUNTERPART'] if v in set(truth.L1) | set(pred)]
        fig.legend(handles=[Line2D([], [], marker='o', color='none', markerfacecolor=colors[v], markersize=5, label=v) for v in labels],
                   loc='outside lower center', ncol=5, fontsize=7, frameon=False)
        files = []
        for suffix in ['png', 'pdf']:
            path = figdir / f'{name}.{suffix}';fig.savefig(path, dpi=180, bbox_inches='tight');files.append(path.name)
        plt.close(fig)
        records.append(dict(condition=condition, terminal_manifest_sha256=sha(source / 'terminal/L00_mean/terminal_manifest.json'),
                            training_executed=tm['training_executed'], dl_status=tm['dl_status'], files=files))
    write_json(figdir / 'manifest.json', dict(status='completed', sample=sample, budget=budget, space=space,
        every_candidate_plotted=True, conditions=records, n_cells=len(truth), display=str(display), display_sha256=sha(display),
        all_cell_denominator=True, source_sha256=sha(__file__), job=os.environ['SLURM_JOB_ID'], completed_at=utc(),
        files={p.name: sha(p) for p in figdir.iterdir() if p.suffix in ['png', '.png', '.pdf']}))
    complete(figdir)
    print('EVALUATION_FIGURES_COMPLETE', sample, budget, space, len(records), flush=True)


if __name__ == '__main__':
    run(sys.argv[1])
