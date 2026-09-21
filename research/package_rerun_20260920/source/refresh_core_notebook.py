#!/usr/bin/env python3
"""Refresh existing concise notebook only after all released-package core fits pass.

Run in SLURM with the campaign's plotting environment. This preserves appended
reviewer sections and changes only core figures/tables and their provenance.
"""
import argparse
import base64
from datetime import datetime, timezone
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CAMPAIGN = ROOT / 'results/hvg_ptc_20260916_v1/package_reference_rerun_20260920'
NOTEBOOK = ROOT / 'notebooks/dgscrna_results.ipynb'
PLOT_SOURCE = ROOT / 'handoff/lfine_compact_20260920/plot_shared_umap.py'


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def require(value, message):
    if not value:
        raise RuntimeError(message)


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def markdown_table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
        '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
        ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gate', required=True, type=Path)
    parser.add_argument('--candidate-only', action='store_true')
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID'), 'Notebook assembly requires SLURM')
    require(not sys.flags.optimize, 'Do not run with Python -O')
    sys.path.insert(0, str(HERE))
    import manage_core
    gate, tasks, root = manage_core.gate_metadata(args.gate.resolve())
    require(root == CAMPAIGN, 'Unexpected campaign root')
    gate_hash = sha(args.gate)
    complete_path = root / 'GBM_CORE_COMPLETE'
    require(complete_path.is_file(), 'Full GBM core acceptance is still pending')
    complete = load(complete_path)
    require(complete['gate_sha256'] == gate_hash and complete['accepted_tasks'] == 726
        and complete['terminal_conditions'] == 139392
        and complete['lfine_valid_threshold_rows'] == 278784, 'Incomplete core acceptance')
    for task in tasks:
        receipt = manage_core.acceptance(root, task, gate_hash, gate['wheel']['sha256'])
        require(receipt == complete['acceptances'][str(task['index'])], 'Core acceptance changed')

    out = root / 'report'
    out.mkdir(parents=True, exist_ok=True)
    evaluation = root / 'evaluation'
    subprocess.run([gate['runtime']['python'], '-s', str(HERE / 'evaluate_core_lfine.py'),
        'aggregate', '--out', str(evaluation)], check=True, cwd='/tmp')
    em = load(evaluation / 'manifest.json')
    require((evaluation / 'COMPLETE').read_text().strip() == sha(evaluation / 'manifest.json'),
        'Evaluation completion checksum mismatch')
    require(em['status'] == 'completed' and em['n_units_evaluated'] == 726
        and em['n_valid_threshold_rows'] == 278784, 'Incomplete cohort evaluation')
    for name, digest in em['outputs'].items():
        require(sha(evaluation / name) == digest, 'Evaluation output changed: ' + name)

    import numpy as np
    import pandas as pd
    import nbformat
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    markers = pd.read_csv(evaluation / 'summary_fullgene_markers.csv')
    hvg = pd.read_csv(evaluation / 'summary_hvg24.csv')
    markers = markers[markers.stage.eq('terminal090')]
    hvg = hvg[hvg.stage.eq('terminal090')]
    require(len(markers) == 64 and len(hvg) == 96, 'Wrong report endpoint roster')
    require(markers.n_samples_unavailable.eq(0).all() and hvg.n_samples_unavailable.eq(0).all(),
        'Unavailable results cannot enter the complete report')
    primary = markers[markers.scope.eq('primary97')].set_index('library')
    worked = markers[markers.scope.eq('TKU3186')].set_index('library')
    overlap = {'CARE_TME', 'BrainAtlas112', 'UNION_all'}
    table = markdown_table(['Marker library', 'TKU3186 F1', 'Primary cohort F1',
        'Coverage', 'DL trained / 97'], [[library + (' †' if library in overlap else ''),
        f'{worked.loc[library, "lfine_macroF1_patient_mean"]:.3f}',
        f'{row.lfine_macroF1_patient_mean:.3f}', f'{row.coverage_patient_mean:.1%}',
        int(row.n_training_executed)] for library, row in primary.iterrows()])
    table += ('\n\n† Source overlaps author-label construction; this is a dependent reference '
        'comparison. Unspecified overlap does not establish independence. Valid no-training '
        'terminal states remain included in F1.')

    figures = out / 'figures'
    figures.mkdir(exist_ok=True)
    plotter = module('packaged_main_figure', PLOT_SOURCE)
    plotter.OUT = figures
    plotter.R_PREP = root / 'core/TKU3186/all/GBM/TKU3186/all'
    plotter.R_ROUTE = plotter.R_PREP / 'UMAP2_HDBSCAN_R'
    truth_path = plotter.NATIVE / 'evaluation_inputs/TKU3186/truth.csv.gz'
    truth_input = load(plotter.NATIVE / 'inputs/TKU3186/input_manifest.json')
    require(sha(truth_path) == truth_input['evaluation_files']['truth.csv.gz'],
        'Author truth used in the main figure changed')
    plotter.main()
    plot_manifest_path = figures / 'native_R_shared_umap_manifest.json'
    plot_manifest = load(plot_manifest_path)
    require(plot_manifest['display_coordinates'] ==
        str((plotter.R_PREP / 'UMAP2.csv').relative_to(ROOT)), 'Main figure used another run')
    for relative_path, digest in plot_manifest['sources'].items():
        require(sha(ROOT / relative_path) == digest, 'Main-figure source changed: ' + relative_path)
    budgets = ['all', 'hvg500', 'hvg1000', 'hvg2000', 'hvg3000', 'hvg5000']
    routes = ['PCA30_SNN', 'PCA30_HDBSCAN_R', 'UMAP2_SNN', 'UMAP2_HDBSCAN_R']
    panels = [hvg[hvg.scope.eq(scope)].pivot(index='budget', columns='route',
        values='lfine_macroF1_patient_mean').loc[budgets, routes]
        for scope in ['TKU3186', 'primary97']]
    require(all(panel.shape == (6, 4) and np.isfinite(panel.values).all() for panel in panels),
        'Incomplete feature/route figure')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.9), constrained_layout=True)
    vmax = max(float(panel.values.max()) for panel in panels)
    for ax, panel, title in zip(axes, panels, ['TKU3186', '97 samples / 55 patients']):
        im = ax.imshow(panel.values, cmap='Blues', vmin=0, vmax=vmax, aspect='auto')
        ax.set_xticks(range(4), ['PCA30\nSNN', 'PCA30\nHDBSCAN', 'UMAP2\nSNN', 'UMAP2\nHDBSCAN'])
        ax.set_yticks(range(6), ['All retained genes', 'HVG 500', 'HVG 1,000',
            'HVG 2,000', 'HVG 3,000', 'HVG 5,000'])
        ax.set_title(title)
        for row in range(6):
            for col in range(4):
                value = panel.iloc[row, col]
                ax.text(col, row, f'{value:.3f}', ha='center', va='center',
                    color='white' if value > .62 * vmax else '#20252b')
        ax.axhline(.5, color='#444444', linewidth=1)
    fig.suptitle('Packaged R + Python DL · Lfine macro-F1 · final 0.90')
    fig.colorbar(im, ax=axes, fraction=.025, pad=.02, label='Macro-F1')
    heatmap = figures / 'native_R_Lfine_hvg_routes.png'
    fig.savefig(heatmap, dpi=170)
    fig.savefig(heatmap.with_suffix('.pdf'))
    plt.close(fig)

    lock_path = ROOT / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/notebook/update.lock'
    with lock_path.open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        before = sha(NOTEBOOK)
        report = nbformat.read(NOTEBOOK, as_version=4)
        require(report.metadata.get('compact_reference_report'), 'Expected existing concise notebook')
        old_sources = [cell.source for cell in report.cells]
        replacements = {}
        def replace(prefix, value):
            positions = [i for i, cell in enumerate(report.cells) if cell.source.startswith(prefix)]
            require(len(positions) == 1, 'Missing/duplicate report section: ' + prefix)
            report.cells[positions[0]].source = value
            replacements[positions[0]] = prefix
        replace('| Marker library |', table)
        replace('# 6 · Robustness', '# 6 · Robustness\n\n'
            'The published package was fitted afresh on all 121 samples × six feature budgets. '
            'All 139,392 terminal conditions passed independent reference checks; both confidence '
            'endpoints were evaluated (278,784 Lfine rows). The fixed primary cohort remains '
            '97 samples / 55 patients; the 24-sample sensitivity subset is not an independent '
            'patient cohort. Reviewer sections below retain their existing results pending '
            'separate package-rerun verification.')
        p = panels[1]
        replace('# 7 · Interpretation', '# 7 · Interpretation\n\n'
            f'With fixed glioma markers and mean cutoff, UMAP/HDBSCAN Lfine macro-F1 is '
            f'**{p.loc["all", "UMAP2_HDBSCAN_R"]:.3f}** with all retained genes and '
            f'**{p.loc["hvg2000", "UMAP2_HDBSCAN_R"]:.3f}** with HVG2000. '
            f'HVG2000/UMAP/SNN scores **{p.loc["hvg2000", "UMAP2_SNN"]:.3f}**. '
            'These are feature/route comparisons; they do not establish uniform optimality.')
        relative = '../results/hvg_ptc_20260916_v1/package_reference_rerun_20260920/'
        replace('# 8 · Data and code', '# 8 · Data and code\n\n'
            'https://github.com/Yiminnn/DGscRNA/releases/tag/v2.0.0rc1\n\n' +
            '\n'.join(f'- [{name}]({relative}{path})' for name, path in [
                ('Complete core metrics', 'evaluation/metrics_all_conditions.csv.gz'),
                ('Marker comparison', 'evaluation/summary_fullgene_markers.csv'),
                ('Feature/route comparison', 'evaluation/summary_hvg24.csv'),
                ('Package acceptance', 'GBM_CORE_COMPLETE'),
                ('Evaluation provenance', 'evaluation/manifest.json')]))
        for name in ['TKU3186_native_R_allgenes_Lfine_marker_contexts.png', heatmap.name]:
            matches = [i for i, cell in enumerate(report.cells) if name in cell.get('attachments', {})]
            require(len(matches) == 1, 'Missing/duplicate core figure: ' + name)
            index = matches[0]
            report.cells[index].attachments[name] = {'image/png': base64.b64encode((figures / name).read_bytes()).decode()}
        # Preserve every other section, especially reviewer additions and PTC.
        require(all(cell.source == old_sources[i] for i, cell in enumerate(report.cells)
            if i not in replacements), 'An unrelated report section changed')
        text = '\n'.join(cell.source for cell in report.cells)
        require('L1' not in text and 'primary hypothesis' not in text.lower(), 'Presentation policy violated')
        report.metadata['packaged_core_rerun'] = dict(gate_sha256=gate_hash,
            release='v2.0.0rc1', core_tasks=726, terminal_conditions=139392,
            lfine_rows=278784, endpoint_displayed='terminal090',
            evaluation_manifest_sha256=sha(evaluation / 'manifest.json'),
            core_completion_sha256=sha(complete_path), builder_sha256=sha(__file__),
            main_figure_manifest_sha256=sha(plot_manifest_path))
        nbformat.validate(report)
        candidate = out / 'dgscrna_results_candidate.ipynb'
        nbformat.write(report, candidate)
        require(sha(NOTEBOOK) == before, 'Concurrent notebook edit detected')
        if not args.candidate_only:
            backup = out / ('notebook_before_package_' + before + '.ipynb')
            if not backup.exists():
                shutil.copy2(NOTEBOOK, backup)
            temporary = NOTEBOOK.with_name('.dgscrna_results.package.tmp')
            shutil.copyfile(candidate, temporary)
            temporary.replace(NOTEBOOK)
        receipt = dict(status='candidate' if args.candidate_only else 'updated',
            updated_at=datetime.now(timezone.utc).isoformat(), job=os.environ['SLURM_JOB_ID'],
            notebook_before_sha256=before, notebook_candidate_sha256=sha(candidate),
            core_gate_sha256=gate_hash, replaced_sections=replacements,
            remaining_sections_preserved=True, PTC_changed=False, HTML_created=False,
            OneDrive_accessed=False, source_hashes={str(path): sha(path) for path in
                [Path(__file__), PLOT_SOURCE, evaluation / 'manifest.json', complete_path,
                 plot_manifest_path, truth_path]},
            figure_hashes={p.name: sha(p) for p in figures.glob('*.png')})
        (out / 'notebook_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
        print(json.dumps(receipt))


if __name__ == '__main__':
    main()
