#!/usr/bin/env python3
"""Frozen patient-fold A1 selection after every new packaged unit is accepted."""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import importlib.util
import json
import os
import shutil
import sys

CODE = Path(__file__).resolve().parent
ROOT = CODE.parents[2]
CAMPAIGN = ROOT / 'results/hvg_ptc_20260916_v1/package_reference_rerun_20260920'
STAGES = {'terminal070', 'terminal090'}


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def require(value, message):
    if not value:
        raise RuntimeError(message)


def receipt(directory):
    directory = Path(directory)
    require((directory / 'COMPLETE').is_file(), 'Completion pending: ' + str(directory))
    require((directory / 'COMPLETE').read_text().strip() == sha(directory / 'manifest.json'),
        'Completion hash changed: ' + str(directory))
    manifest = load(directory / 'manifest.json')
    for name, digest in manifest['outputs'].items():
        require(sha(directory / name) == digest, 'Output changed: ' + name)
    return manifest


def source_contract():
    contract = load(CODE / 'SOURCE_PROTOCOL_MANIFEST.json')
    for path, digest in contract['source_hashes'].items():
        require(sha(path) == digest, 'Original selection source changed: ' + path)
    require(sha(CODE / 'numerical.py') == contract['adapted_numerical_sha256'],
        'Scientific selection body changed')
    require(sha(CODE / 'canonical_inference.py') == contract['independent_inference']['sha256'],
        'Independent canonical inference verifier changed')
    for name, original in [('lfine.json', 'embedding_lfine_v1/protocol.json'),
            ('embedding.json', 'protocol/embedding.json'),
            ('selection.json', 'protocol/embedding_selection.json'),
            ('patient_folds.csv', 'protocol/embedding_patient_folds.csv')]:
        source = ROOT / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920' / original
        require(sha(CODE / 'protocol' / name) == sha(source), 'Copied protocol changed')
    return contract


def prerequisites(args):
    """Metadata and checksums only: refuse partial campaigns before loading metrics."""
    require(os.environ.get('SLURM_JOB_ID'), 'Scientific selection requires SLURM')
    require(not sys.flags.optimize, 'Do not run with Python -O')
    contract = source_contract()
    complete_path = CAMPAIGN / 'GBM_EXTENSIONS_COMPLETE'
    require(complete_path.is_file(), 'All 1,694 packaged A1 units are not yet accepted')
    completion = load(complete_path)
    require(completion['accepted_tasks'] == 2617 and
        completion['campaign_gate_sha256'] == sha(args.campaign), 'Wrong extension completion')
    core_path = CAMPAIGN / 'GBM_CORE_COMPLETE'
    require(core_path.is_file(), 'Complete new core anchors are required')
    core_completion = load(core_path)
    campaign_gate = load(args.campaign)
    release_gate = campaign_gate['release_gate']
    require(sha(release_gate['path']) == release_gate['sha256'] and
        core_completion['gate_sha256'] == release_gate['sha256'], 'Core anchors belong to another release')
    require(core_completion['accepted_tasks'] == 726 and
        core_completion['lfine_valid_threshold_rows'] == 278784, 'Incomplete core anchors')
    aggregate = receipt(args.aggregate)
    require(aggregate['status'] == 'completed' and aggregate['tasks'] == 2617 and
        aggregate['completion_sha256'] == sha(complete_path) and
        aggregate['campaign_gate_sha256'] == sha(args.campaign), 'Unaccepted extension aggregation')
    core = receipt(CAMPAIGN / 'evaluation')
    require(core['status'] == 'completed' and core['n_units_evaluated'] == 726 and
        core['n_valid_threshold_rows'] == 278784, 'Incomplete new core evaluation')
    require(core['frozen_semantics'] == aggregate['frozen_semantics'], 'Lfine endpoint differs')
    require(not args.out.exists(), 'Preserve previous summary attempts')
    inputs = {str(path): sha(path) for path in [args.campaign, complete_path, core_path,
        args.aggregate / 'manifest.json', CAMPAIGN / 'evaluation/manifest.json',
        args.aggregate / 'metrics_all_extensions.csv.gz',
        CAMPAIGN / 'evaluation/metrics_hvg24.csv.gz', CODE / 'SOURCE_PROTOCOL_MANIFEST.json',
        CODE / 'numerical.py', CODE / 'canonical_inference.py', Path(__file__)]}
    for path, digest in aggregate['input_hashes'].items():
        require(sha(path) == digest, 'Accepted extension input changed: ' + path)
    return contract, inputs


def compact_text(statistics):
    subset = statistics[statistics.cohort.eq('primary97')]
    rows = ['## DGCyTOF-style extension', '',
        'New packaged runs: 121 samples × two HVG budgets × seven spaces × 13 partitions. '
        'The primary cohort totals 97 samples / 55 patients; scores use terminal 0.90 Lfine compatible-target-set F1.', '',
        '| Representation / clusterer | HVG2000 F1 (Δ) | HVG5000 F1 (Δ) |',
        '|---|---:|---:|']
    for space in ['noDR', 'PCA2', 'FA2', 'ICA2', 'Isomap2', 'UMAP2', 'TSNE2']:
        for method in ['KMeans', 'GMM', 'HDBSCAN_R']:
            cells = []
            for budget in ['hvg2000', 'hvg5000']:
                match = subset[subset.space.eq(space) & subset.method.eq(method) & subset.budget.eq(budget)]
                require(len(match) == 1, 'Incomplete compact comparison table')
                row = match.iloc[0]
                cells.append(f'{row.candidate_mean_F1:.3f} ({row.mean_delta:+.3f})')
            label = 'ICA2 adaptive' if space == 'ICA2' else space
            rows.append('| ' + label + ' / ' + method.replace('_R', '') + ' | ' + ' | '.join(cells) + ' |')
    anchor = [float(subset[subset.budget.eq(b)].anchor_mean_F1.iloc[0]) for b in ['hvg2000', 'hvg5000']]
    rows += ['', f'Δ compares the original PCA30 → UMAP → HDBSCAN anchor (F1 {anchor[0]:.3f} / {anchor[1]:.3f}). '
        'The UMAP comparator uses scaled HVGs directly. K is chosen on training patients only '
        '(five fixed folds; exact ties choose smaller K); HDBSCAN is fixed.', '',
        'All 21 contrasts per budget are retained. Patient bootstrap confidence intervals and Holm-adjusted '
        'Wilcoxon tests are conditional on the selected predictions; they do not establish uniform optimality. '
        'Evaluable sample/patient denominators accompany every linked contrast.', '',
        '[All contrasts, confidence intervals and denominators](../results/hvg_ptc_20260916_v1/'
        'package_reference_rerun_20260920/evaluation/embedding_selected/paired_vs_original_anchor.csv) · '
        '[Fold choices](../results/hvg_ptc_20260916_v1/package_reference_rerun_20260920/'
        'evaluation/embedding_selected/patient_fold_K_choices.csv)']
    return '\n'.join(rows)


def update_notebook(out, text, update):
    import nbformat
    notebook = ROOT / 'notebooks/dgscrna_results.ipynb'
    lockpath = ROOT / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/notebook/update.lock'
    with lockpath.open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        before = sha(notebook)
        report = nbformat.read(notebook, as_version=4)
        require(report.metadata.get('compact_reference_report'), 'Expected the existing concise notebook')
        positions = [i for i, cell in enumerate(report.cells)
            if cell.source.startswith('## DGCyTOF-style extension')]
        require(len(positions) == 1, 'Missing/duplicate existing DGCyTOF section')
        untouched = {i: json.dumps(cell, sort_keys=True) for i, cell in enumerate(report.cells) if i != positions[0]}
        report.cells[positions[0]].source = text
        report.metadata['packaged_A1_selection'] = dict(summary_sha256=sha(out / 'manifest.json'),
            terminal_only=True, primary_samples=97, primary_patients=55, selection_rows=420,
            contrasts=84, all_candidate_conditions_retained=True)
        for index, value in untouched.items():
            require(json.dumps(report.cells[index], sort_keys=True) == value, 'Unrelated notebook cell changed')
        nbformat.validate(report)
        candidate = out / 'dgscrna_results_candidate.ipynb'
        nbformat.write(report, candidate)
        require(sha(notebook) == before, 'Concurrent notebook change detected')
        if update:
            shutil.copy2(notebook, out / ('notebook_before_' + before + '.ipynb'))
            temporary = notebook.with_name('.dgscrna_results.A1.tmp')
            shutil.copyfile(candidate, temporary)
            temporary.replace(notebook)
        record = dict(status='updated' if update else 'candidate', notebook_before_sha256=before,
            notebook_candidate_sha256=sha(candidate), replaced_cell=positions[0],
            remaining_cells_preserved=True, PTC_changed=False, HTML_created=False, OneDrive_accessed=False)
        (out / 'notebook_receipt.json').write_text(json.dumps(record, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', type=Path, required=True, help='Unified extension campaign gate')
    parser.add_argument('--aggregate', type=Path, required=True)
    parser.add_argument('--out', type=Path, default=CAMPAIGN / 'evaluation/embedding_selected')
    parser.add_argument('--update-notebook', action='store_true')
    args = parser.parse_args()
    args.campaign, args.aggregate, args.out = args.campaign.resolve(), args.aggregate.resolve(), args.out.resolve()
    require(args.out == CAMPAIGN / 'evaluation/embedding_selected', 'Use the canonical new A1 summary directory')
    contract, inputs = prerequisites(args)
    import pandas as pd
    spec_module = importlib.util.spec_from_file_location('frozen_A1_selection', CODE / 'numerical.py')
    numerical = importlib.util.module_from_spec(spec_module)
    spec_module.loader.exec_module(numerical)
    data = pd.read_csv(args.aggregate / 'metrics_all_extensions.csv.gz', dtype={'sample': str, 'patient': str, 'cutoff': str})
    data = data[data.family.eq('seven_space_control')].copy()
    require(len(data) == 44044 and data.campaign_task.nunique() == 1694, 'Incomplete A1 threshold roster')
    require(data.terminal_valid.all() and data.status.eq('completed').all(), 'Invalid A1 terminal state')
    require(data.library.eq('CM2_glioma_other').all() and data.cutoff.eq('mean').all(), 'Marker context changed')
    spec = load(CODE / 'protocol/embedding.json')
    protocol = load(CODE / 'protocol/lfine.json')
    folds = pd.read_csv(CODE / 'protocol/patient_folds.csv', dtype={'sample': str, 'patient': str})
    require(sha(CODE / 'protocol/patient_folds.csv') == protocol['patient_folds_sha256'], 'Patient folds changed')
    require(len(folds) == 121 and folds['sample'].nunique() == 121 and folds.primary.sum() == 97 and
        folds.groupby('patient').fold.nunique().eq(1).all(), 'Patient/fold roster mismatch')
    data['space'] = data.configuration
    partitions = [(f'{method}_K{k:02d}', method, k) for method in ['KMeans', 'GMM'] for k in spec['K']] + [('HDBSCAN_R', 'HDBSCAN_R', 0)]
    routes = {f'{space}_{name}': (space, method, k) for space in spec['spaces'] for name, method, k in partitions}
    require(set(data.route) == set(routes), 'Candidate route roster changed')
    require(all(routes[route][0] == space for route, space in zip(data.route, data.space)), 'Route/space mismatch')
    data['method'] = data.route.map(lambda route: routes[route][1])
    data['k'] = data.route.map(lambda route: routes[route][2])
    data['lfine_metric_eligible'] = data.lfine_n_classes.gt(0)
    keys = ['sample', 'budget', 'space', 'method', 'k', 'stage']
    require(set(data.stage) == STAGES and not data.duplicated(keys).any(), 'Duplicate/wrong candidate endpoints')
    require(data.groupby(['sample', 'budget', 'space']).size().eq(26).all(), 'Incomplete representation endpoint roster')
    anchors = pd.read_csv(CAMPAIGN / 'evaluation/metrics_hvg24.csv.gz', dtype={'sample': str, 'patient': str, 'cutoff': str})
    anchors = anchors[anchors.budget.isin(spec['budgets']) & anchors.route.eq('UMAP2_HDBSCAN_R')].copy()
    require(len(anchors) == 484 and anchors.terminal_valid.all() and set(anchors.stage) == STAGES,
        'Incomplete same-budget original anchors')
    require(anchors.library.eq('CM2_glioma_other').all() and anchors.cutoff.eq('mean').all(), 'Anchor context changed')
    anchors['space'], anchors['method'], anchors['k'] = 'PCA30_UMAP2', 'HDBSCAN_R', 0
    anchors['lfine_metric_eligible'] = anchors.lfine_n_classes.gt(0)
    require(data.lfine_macroF1.notna().equals(data.lfine_metric_eligible) and
        anchors.lfine_macroF1.notna().equals(anchors.lfine_metric_eligible), 'Candidate-dependent metric eligibility')
    eligibility = {}
    fold_lookup = folds.set_index('sample')
    for sample, group in data.groupby('sample'):
        metadata = group[['patient', 'primary', 'n_cells', 'lfine_n_classes', 'truth_sha256']].drop_duplicates()
        require(len(metadata) == 1, 'Candidate-dependent truth/denominator')
        row = metadata.iloc[0]
        require(row.patient == fold_lookup.loc[sample, 'patient'] and bool(row.primary) == bool(fold_lookup.loc[sample, 'primary']), 'Patient cohort changed')
        reference = anchors[anchors['sample'].eq(sample)]
        require(reference.truth_sha256.eq(row.truth_sha256).all() and
            reference.lfine_n_classes.eq(row.lfine_n_classes).all(), 'Candidate/anchor truth mismatch')
        eligibility[sample] = dict(sample=sample, patient=row.patient, primary=bool(row.primary),
            n_cells=int(row.n_cells), n_scored_classes=int(row.lfine_n_classes),
            lfine_metric_eligible=bool(row.lfine_n_classes), truth_sha256=row.truth_sha256)
    require(len(eligibility) == 121 and data.patient.nunique() == 59, 'Incomplete sample/patient roster')
    choices, selected, patients, stats = numerical.select(data, anchors, eligibility, folds, spec, protocol)
    independent = numerical.independently_select(data, anchors, eligibility, folds, spec, protocol)
    require(len(choices) == 420 and len(stats) == 84, 'Incomplete selected comparisons')
    args.out.mkdir(parents=True)
    inference_spec = importlib.util.spec_from_file_location('A1_canonical_inference', CODE / 'canonical_inference.py')
    inference = importlib.util.module_from_spec(inference_spec)
    inference_spec.loader.exec_module(inference)
    inference_proof = inference.verify(numerical, (choices, selected, patients, stats), independent,
        anchors, spec, args.out / 'inference_verification')
    for path, digest in inputs.items():
        require(sha(path) == digest, 'Input changed during selection: ' + path)
    source_contract()
    tables = dict(all_candidate_metrics=data, anchor_metrics=anchors, selected_sample_metrics=selected,
        patient_metrics=patients, patient_fold_K_choices=choices, paired_vs_original_anchor=stats)
    for name, frame in tables.items():
        suffix = '.csv' if name in {'patient_fold_K_choices', 'paired_vs_original_anchor'} else '.csv.gz'
        frame.to_csv(args.out / (name + suffix), index=False)
    (args.out / 'sample_eligibility.json').write_text(json.dumps(eligibility, indent=2) + '\n')
    text = compact_text(stats)
    (args.out / 'notebook_section.md').write_text(text + '\n')
    manifest = dict(status='passed_independent_selection', n_representations=1694,
        n_candidate_rows=44044, n_anchor_rows=484, n_selection_rows=420, n_paired_contrasts=84,
        terminal_only=True, fitting_performed=False, all_candidates_retained=True,
        selection='Frozen training-patient terminal090 Lfine selection, exact tie smaller K; held-out patients never select K',
        inference='Retrospective patient-label-held-out K choice; inference conditional on selected predictions',
        independent_selection_and_statistics_passed=True, source_contract=contract,
        canonical_inference_verification=inference_proof,
        source_sha256=sha(__file__), inputs=inputs,
        outputs={str(p.relative_to(args.out)): sha(p) for p in args.out.rglob('*') if p.is_file()},
        job=os.environ['SLURM_JOB_ID'], completed_at=datetime.now(timezone.utc).isoformat())
    (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=2, allow_nan=False) + '\n')
    (args.out / 'COMPLETE').write_text(sha(args.out / 'manifest.json') + '\n')
    update_notebook(args.out, text, args.update_notebook)
    print(json.dumps({key: manifest[key] for key in ['status', 'n_representations', 'n_selection_rows', 'n_paired_contrasts']}))


if __name__ == '__main__':
    main()
