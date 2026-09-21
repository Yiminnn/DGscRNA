"""Evaluate only frozen predictions, and plot every candidate on fixed coordinates."""
from pathlib import Path
import json
import os
import shutil
import sys
CODE = Path('/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/reviewer_completion_20260920/source_snapshots/embedding_v8')
sys.path.insert(0, str(CODE))
sys.path.insert(0, str(CODE / 'legacy'))
import run as core
from common import OUT as REFERENCE, L1, require_slurm, write_json, sha, checked, complete, utc


def render(directory):
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
    core.cache_compatibility.own_source(core)
    cfg = json.loads((directory / 'config.json').read_text())
    sample, budget, space = cfg['sample'], cfg['budget'], cfg['space']
    assert space=='ICA2', 'This versioned renderer only adds adaptive ICA figures'
    policy=CODE.parent.parent/'protocol/embedding_convergence_repair_20260920_v2.json'
    policy_data=json.loads(policy.read_text())
    assert policy_data['policy_id']=='ICA2_parallel_caps_then_deflation_v2'
    representation=json.loads((directory/'representation.json').read_text())
    solver=representation.get('actual_solver',representation['params']['algorithm'])
    cap=representation.get('actual_iteration_cap',representation['params']['max_iter'])
    space_display='ICA2 adaptive'
    figdir=directory/'figures_adaptive_v1'
    native_cells = Path(cfg['prep']) / 'cells.csv'
    assert (directory/'cells.csv').exists()
    assert sha(directory / 'cells.csv') == sha(native_cells)
    assert checked(directory/'evaluation')
    if checked(figdir):
        verify_manifest(directory)
        return
    previous=json.loads((directory/'fit_manifest.json').read_text())
    core.cache_compatibility.validate_completed_fit(core,cfg,previous)
    evaluation_manifest=json.loads((directory/'evaluation/manifest.json').read_text())
    assert sha(directory/'evaluation/metrics.csv')==evaluation_manifest['outputs']['metrics.csv']
    metrics = pd.read_csv(directory / 'evaluation/metrics.csv', dtype={'cutoff': str})
    assert len(metrics) == 3 * len(cfg['conditions'])
    display = REFERENCE / 'GBM' / sample / 'hvg2000/UMAP2.csv'
    assert sha(REFERENCE/'evaluation_inputs'/sample/'truth.csv.gz')==evaluation_manifest['truth_sha256']
    assert sha(REFERENCE/'markers/panel_L1_mapping.csv')==evaluation_manifest['mapping_sha256']
    frame = pd.read_csv(display, index_col=0)
    coords = frame.to_numpy()
    truth = pd.read_csv(REFERENCE / 'evaluation_inputs' / sample / 'truth.csv.gz', dtype=str, keep_default_na=False)
    assert list(frame.index) == list(truth.cell_id)
    assert list(pd.read_csv(directory/'cells.csv',dtype=str,keep_default_na=False).cell_id)==list(truth.cell_id)
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

    figdir.mkdir(exist_ok=True)
    records = []
    for condition in cfg['conditions']:
        source = Path(condition['dest']);name = source.name
        assert sha(source/'terminal/L00_mean/terminal_manifest.json')==evaluation_manifest['terminal_manifests'][name+'/terminal/L00_mean']
        tm = json.loads((source / 'terminal/L00_mean/terminal_manifest.json').read_text())
        cl = pd.read_csv(source / 'clusters.csv', dtype=str, keep_default_na=False)
        pp = pd.read_csv(source / 'terminal/L00_mean/predictions.csv.gz', dtype=str, keep_default_na=False)
        assert list(cl.cell_id) == list(pp.cell_id) == list(truth.cell_id)
        pred = [lookup.get(v, 'Unknown' if v in ['Unknown', 'Undecided', 'Noise', ''] else 'UNMAPPABLE') for v in pp.final090]
        row = metrics[(metrics.route == name) & (metrics.stage == 'terminal090')].iloc[0]
        fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.3), layout='constrained')
        panel(axes[0], truth.L1, 'Original author L1')
        panel(axes[1], cl.cluster, f'{space_display} / {name}: partition', True, condition['method'] == 'HDBSCAN_R')
        stage_title = 'Final DL threshold 0.90' if tm['training_executed'] else 'Terminal result: DL not executed'
        panel(axes[2], pred, f'{stage_title}\n{tm["dl_status"]}\nmacro-F1 {row.macroF1_present:.3f}; coverage {row.coverage:.1%}')
        fig.suptitle(f'{sample} | {budget} | {space_display} / {name}\nParallel, then deflation fallback; actual {solver}, cap {cap}\nCM2_glioma_other / mean; shared display coordinates', fontsize=10)
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
        space_display=space_display,comparator_policy=policy_data['comparator_name'],
        policy_sha256=sha(policy),actual_solver=solver,actual_iteration_cap=cap,
        representation_manifest_sha256=sha(directory/'representation.json'),
        evaluation_manifest_sha256=sha(directory/'evaluation/manifest.json'),
        fit_manifest_sha256=sha(directory/'fit_manifest.json'),
        presentation_source_manifest_sha256=sha(Path(__file__).resolve().parent/'SOURCE_MANIFEST.json'),
        fit_consumer_source_manifest_sha256=sha(CODE/'SOURCE_MANIFEST.json'),
        canonical_parallel5000_status=representation.get('canonical_parallel5000_status','converged'),
        fallback_used=representation.get('fallback_used',False),
        old_figures_preserved=True,metrics_recomputed=False,
        every_candidate_plotted=True, conditions=records, n_cells=len(truth), display=str(display), display_sha256=sha(display),
        all_cell_denominator=True, source_sha256=sha(__file__), job=os.environ['SLURM_JOB_ID'], completed_at=utc(),
        files={p.name: sha(p) for p in figdir.iterdir() if p.suffix in ['png', '.png', '.pdf']}))
    complete(figdir)
    print('EVALUATION_FIGURES_COMPLETE', sample, budget, space, len(records), flush=True)


def own_source():
    manifest=json.loads((Path(__file__).resolve().parent/'SOURCE_MANIFEST.json').read_text())
    for name,digest in manifest.items():
        assert sha(Path(__file__).resolve().parent/name)==digest, 'Changed presentation source'


def verify_manifest(directory):
    directory=Path(directory);figdir=directory/'figures_adaptive_v1'
    assert checked(figdir)
    m=json.loads((figdir/'manifest.json').read_text())
    for field,path in [('source_sha256',Path(__file__)),
                       ('presentation_source_manifest_sha256',Path(__file__).resolve().parent/'SOURCE_MANIFEST.json'),
                       ('representation_manifest_sha256',directory/'representation.json'),
                       ('evaluation_manifest_sha256',directory/'evaluation/manifest.json'),
                       ('fit_manifest_sha256',directory/'fit_manifest.json'),
                       ('policy_sha256',CODE.parent.parent/'protocol/embedding_convergence_repair_20260920_v2.json')]:
        assert m[field]==sha(path), f'Changed adaptive figure input/source: {field}'
    assert m['space_display']=='ICA2 adaptive' and m['every_candidate_plotted'] is True
    assert len(m['conditions'])==13 and len(m['files'])==26
    expected={Path(c['condition']['dest']).name+'.'+ext for c in m['conditions'] for ext in ('png','pdf')}
    assert set(m['files'])==expected
    for name,digest in m['files'].items():assert sha(figdir/name)==digest
    return m


def run(directory):
    require_slurm();own_source()
    directory=Path(directory)
    cfg=json.loads((directory/'config.json').read_text())
    if cfg['space']!='ICA2':
        print('ADAPTIVE_PRESENTATION_NOT_APPLICABLE',cfg['sample'],cfg['budget'],cfg['space'],flush=True)
        return
    lock=directory/'figures_adaptive_v1.render.lock'
    lock.mkdir()  # Atomic across nodes; never guess that a live/stale owner is gone.
    write_json(lock/'owner.json',dict(job=os.environ['SLURM_JOB_ID'],pid=os.getpid(),source_sha256=sha(__file__)))
    try:
        render(directory)
    finally:
        (lock/'owner.json').unlink();lock.rmdir()


def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('directory',nargs='?')
    parser.add_argument('--tasks')
    args=parser.parse_args()
    assert bool(args.directory)!=bool(args.tasks), 'Pass directory OR --tasks tasks.json'
    if args.tasks:
        tasks=json.loads(Path(args.tasks).read_text())
        task=tasks[int(os.environ.get('SLURM_ARRAY_TASK_ID','0'))]
        for space in task.get('spaces',[task.get('space')]):
            if space!='ICA2':
                print('ADAPTIVE_PRESENTATION_NOT_APPLICABLE',task['sample'],task['budget'],space,flush=True)
                continue
            run(core.OUTPUT/task['sample']/task['budget']/space)
    else:run(args.directory)


if __name__=='__main__':main()
