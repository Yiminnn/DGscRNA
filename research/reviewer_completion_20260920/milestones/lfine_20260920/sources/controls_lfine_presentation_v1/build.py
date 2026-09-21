"""Concise B artifacts from independently accepted terminal Lfine results."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os

assert os.environ.get('SLURM_JOB_ID'), 'Scientific summaries and plots require SLURM'
ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
CAMP = ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
SRC = CAMP/'controls_lfine_v1'
OUT = CAMP/'controls_lfine_presentation_v1'


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


proof = read(SRC/'independent_verification/validation.json')
full = read(SRC/'full/validation.json')
assert proof['status'] == 'passed' and proof['n_rows'] == 408
assert proof['full_validation_sha256'] == sha(SRC/'full/validation.json')
assert full['status'] == 'passed' and full['no_best_selection']
for name, digest in full['files'].items():
    assert sha(SRC/'full'/name) == digest, name
inputs = {str(path.relative_to(ROOT)): sha(path) for path in [
    SRC/'independent_verification/validation.json', SRC/'full/validation.json',
    SRC/'full/terminal_metrics.csv.gz', CAMP/'controls/tasks.json']}
assert not OUT.exists(), 'Preserve prior presentation; create a new version to revise it'
OUT.mkdir()

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

frame = pd.read_csv(SRC/'full/terminal_metrics.csv.gz')
assert len(frame) == 408 and frame.row_key.is_unique and frame.terminal_valid.all()
samples = ['TKU4163', 'NL022', 'SN040']
budgets = ['hvg2000', 'hvg5000']
libraries = ['CM2_glioma_other', 'CM2_primary_all_context']
colors = dict(zip(samples, ['#0072B2', '#D55E00', '#009E73']))
styles = dict(zip(libraries, ['-', '--']))
sweeps = [
    ('PCA30_SNN', 'snn_k', 20, 'PCA–SNN · SNN k', [10, 20, 40]),
    ('UMAP2_SNN', 'snn_k', 20, 'UMAP–SNN · SNN k', [10, 20, 40]),
    ('UMAP2_SNN', 'umap_neighbors', 30, 'UMAP–SNN · UMAP neighbors', [15, 30, 60]),
    ('UMAP2_HDBSCAN_R', 'umap_neighbors', 30, 'UMAP–HDBSCAN · UMAP neighbors', [15, 30, 60]),
]
neighbors = frame[(frame.task == 'neighbors') & (frame.stage == 'terminal090')].copy()
assert len(neighbors) == 132
plt.rcParams.update({'font.size': 8, 'font.family': 'DejaVu Sans',
                     'pdf.fonttype': 42, 'axes.spines.top': False,
                     'axes.spines.right': False})
fig, axes = plt.subplots(2, 4, figsize=(11.8, 5.7), sharey=True)
points, descriptions = [], []
for i, budget in enumerate(budgets):
    for j, (route, parameter, default, title, ticks) in enumerate(sweeps):
        ax = axes[i, j]
        other = 'umap_neighbors' if parameter == 'snn_k' else 'snn_k'
        other_value = 30 if other == 'umap_neighbors' else 20
        ax.axvline(default, color='#bbbbbb', linewidth=.8, zorder=0)
        for sample in samples:
            for library in libraries:
                rows = neighbors[(neighbors['sample'] == sample) &
                                 (neighbors.budget == budget) &
                                 (neighbors.route == route) &
                                 (neighbors.library == library) &
                                 (neighbors[other] == other_value)].sort_values(parameter)
                assert len(rows) == 3 and rows[parameter].tolist() == ticks
                x, y = rows[parameter].to_numpy(), rows.lfine_macroF1.to_numpy()
                assert np.isfinite(y).all()
                line, = ax.plot(x, y, marker='o', markersize=3.5, linewidth=1.25,
                                color=colors[sample], linestyle=styles[library])
                assert np.array_equal(line.get_xdata(), x) and np.array_equal(line.get_ydata(), y)
                plotted = rows.copy()
                plotted['panel'] = f'{budget}/{route}/{parameter}'
                plotted['x'] = x
                plotted['parameter'] = parameter
                points.append(plotted)
                baseline = rows[rows[parameter] == default].iloc[0]
                maximum = float(rows.lfine_macroF1.max())
                relation = 'lower_than_an_alternative' if maximum > baseline.lfine_macroF1 + 1e-12 else (
                    'tied_maximum' if (np.abs(rows.lfine_macroF1 - maximum) <= 1e-12).sum() > 1 else 'unique_maximum')
                descriptions.append(dict(sample=sample, budget=budget, route=route,
                    library=library, parameter=parameter, default=default,
                    default_lfine_macroF1=baseline.lfine_macroF1,
                    default_coverage=baseline.coverage,
                    tested_min_lfine_macroF1=float(rows.lfine_macroF1.min()),
                    tested_max_lfine_macroF1=maximum,
                    default_relation=relation,
                    values=json.dumps([dict(setting=int(r[parameter]), lfine_macroF1=float(r.lfine_macroF1),
                        coverage=float(r.coverage), row_key=r.row_key) for _, r in rows.iterrows()]),
                    interpretation='Descriptive within this three-setting pilot sweep; no selected production parameter'))
        ax.set_title(title, fontsize=8.5)
        ax.set_xticks(ticks)
        ax.set_ylim(-.025, 1.025)
        ax.set_yticks([0, .25, .5, .75, 1])
        ax.grid(axis='y', color='#dddddd', linewidth=.45)
        ax.set_axisbelow(True)
        if j == 0:
            ax.set_ylabel(budget.upper()+'\nLfine compatibility F1')
fig.suptitle('GBM neighborhood sensitivity', fontsize=12, y=.98)
fig.text(.5, .926, 'Original R → terminal DL · fixed marker sets · confidence threshold 0.90',
         ha='center', fontsize=9)
handles = [Line2D([0], [0], color=colors[s], marker='o', label=s) for s in samples]
handles += [Line2D([0], [0], color='#444444', linestyle=styles[l],
                  label={'CM2_glioma_other': 'Glioma/other markers',
                         'CM2_primary_all_context': 'All-context markers'}[l]) for l in libraries]
handles += [Line2D([0], [0], color='#bbbbbb', label='Source default (vertical line)')]
fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(.5, .025))
fig.subplots_adjust(top=.86, bottom=.19, left=.065, right=.985, hspace=.38, wspace=.16)
for ext in ['png', 'pdf']:
    fig.savefig(OUT/f'neighbor_lfine.{ext}', dpi=300)
plt.close(fig)
plot_points = pd.concat(points, ignore_index=True)
assert len(plot_points) == 144 and plot_points.row_key.nunique() == 132
assert set(plot_points.row_key) == set(neighbors.row_key)
plot_points.to_csv(OUT/'neighbor_plot_points.csv.gz', index=False)
pd.DataFrame(descriptions).to_csv(OUT/'neighbor_descriptive_summary.csv', index=False)

learning = frame[frame.task == 'learning'].copy()
assert len(learning) == 144
group = learning.groupby(['sample', 'budget', 'stage', 'epochs'], sort=False)
assert group.size().eq(3).all()
assert group.model_seed.apply(lambda s: set(s) == {0, 1, 42}).all()
checkpoints = group.agg(lfine_macroF1_mean=('lfine_macroF1', 'mean'),
    lfine_macroF1_seed_SD=('lfine_macroF1', 'std'), coverage_mean=('coverage', 'mean'),
    coverage_seed_SD=('coverage', 'std'), known_training_classes=('n_training_classes', 'min'),
    n_seeds=('model_seed', 'nunique')).reset_index()
checkpoints.to_csv(OUT/'checkpoint_lfine_all_thresholds.csv', index=False)
primary = checkpoints[checkpoints.stage == 'terminal090'].copy()
assert len(primary) == 24
lines = ['# GBM parameter sensitivity — Lfine', '',
    'Three size-selected GBM samples; two HVG budgets. These are descriptive pilot controls.', '',
    'The figure includes both fixed marker contexts and every unique neighbor condition at threshold 0.90. '
    'The source-default UMAP–SNN setting is displayed in both applicable one-parameter sweeps; it is one saved result. '
    'Scores use compatible fine-label target sets, not one-to-one fine-subtype predictions.', '',
    '## Training checkpoints', '',
    'Glioma/other markers; terminal 0.90. Each entry is mean F1 / called coverage across initialization seeds 0, 1, 42. '
    'All checkpoints remain reported; no biological-truth early stopping or outcome-based checkpoint replacement.', '',
    '| Sample | HVG | Epoch 5 | Epoch 10 (reference) | Epoch 20 | Epoch 30 | Known seed classes |',
    '|---|---:|---:|---:|---:|---:|---:|']
for sample in samples:
    for budget in budgets:
        part = primary[(primary['sample'] == sample) & (primary.budget == budget)].set_index('epochs')
        assert set(part.index) == {5, 10, 20, 30}
        values = [f'{part.loc[e,"lfine_macroF1_mean"]:.3f} / {part.loc[e,"coverage_mean"]:.3f}' for e in [5, 10, 20, 30]]
        classes = int(part.known_training_classes.min())
        lines.append('| '+' | '.join([sample, budget[3:], *values, str(classes)])+' |')
lines += ['', 'Full seed SD, threshold 0.70 and 0.90 results, and descriptive default-versus-grid ranges are in the linked CSVs. '
    'The TKU4163/HVG5000 learning fits have a single known seed class; their training curves are not multiclass biological validation. '
    'Pseudo-label validation-loss curves remain in the original controls archive. '
    'No p-value, population-wide optimum, or all-gene fit is inferred from these pilots.', '']
(OUT/'REPORT.md').write_text('\n'.join(lines))
counts = pd.DataFrame(descriptions).groupby(['library', 'default_relation']).size()
summary = {'status': 'completed', 'created_at': datetime.now(timezone.utc).isoformat(),
    'job': os.environ['SLURM_JOB_ID'], 'step': os.environ.get('SLURM_STEP_ID'),
    'source_sha256': sha(__file__), 'inputs': inputs,
    'n_terminal_input_rows': 408, 'n_neighbor_unique_displayed_rows': 132,
    'n_neighbor_displayed_points': 144, 'n_descriptive_sweep_families': 48,
    'n_checkpoint_aggregate_rows': len(checkpoints), 'model_fits': 0,
    'all_neighbor_conditions_in_figure': True, 'canonical_notebook_written': False,
    'endpoint': 'Lfine compatible-target-set metrics',
    'descriptive_default_relations': [{'library': k[0], 'relation': k[1], 'n': int(v)} for k, v in counts.items()],
    'not_independent_replicates': 'Sweeps reuse three size-selected samples and sometimes one default point',
    'files': {p.name: sha(p) for p in OUT.iterdir() if p.is_file()}}
for name, digest in inputs.items():
    assert sha(ROOT/name) == digest, name
(OUT/'manifest.json').write_text(json.dumps(summary, indent=2)+'\n')
(OUT/'COMPLETE').write_text(sha(OUT/'manifest.json')+'\n')
print(json.dumps(summary, indent=2), flush=True)
