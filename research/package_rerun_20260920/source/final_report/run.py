#!/usr/bin/env python3
"""Refresh the existing B/A2/C report only after complete packaged GBM acceptance."""
from pathlib import Path
import argparse
import base64
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
import shutil
import sys

CODE = Path(__file__).resolve().parent
ROOT = CODE.parents[2]
CAMP = ROOT / 'results/hvg_ptc_20260916_v1/package_reference_rerun_20260920'
OLD = ROOT / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
STAGES = {'terminal070', 'terminal090'}


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def need(value, message):
    if not value:
        raise RuntimeError(message)


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def checked(directory, inputs, status=None):
    directory = Path(directory)
    path = directory / 'manifest.json'
    need((directory / 'COMPLETE').read_text().strip() == sha(path), 'Completion checksum differs: ' + str(directory))
    value = load(path)
    if status is not None:
        need(value['status'] == status, 'Incomplete source: ' + str(directory))
    inputs[str(path)] = sha(path)
    for name, digest in value.get('outputs', value.get('files', {})).items():
        source = directory / name
        need(sha(source) == digest, 'Changed completed source: ' + str(source))
        inputs[str(source)] = digest
    return value


def prerequisites(args):
    """No pandas, numerical source imports, scientific reads or writes before gates."""
    need(os.environ.get('SLURM_JOB_ID'), 'Report revalidation requires SLURM')
    need(not sys.flags.optimize, 'Do not run with Python -O')
    extension_path, core_path = CAMP / 'GBM_EXTENSIONS_COMPLETE', CAMP / 'GBM_CORE_COMPLETE'
    need(extension_path.is_file(), 'Full GBM extension completion is pending')
    need(core_path.is_file(), 'Full GBM core completion is pending')
    gate = load(args.campaign)
    need(gate['scope'] == 'full_2617' and gate['task_count'] == 2617 and Path(gate['output_root']) == CAMP,
         'Unexpected campaign scope')
    need(sha(gate['release_gate']['path']) == gate['release_gate']['sha256'], 'Changed core release gate')
    ext, core = load(extension_path), load(core_path)
    need(ext['accepted_tasks'] == 2617 and ext['campaign_gate_sha256'] == sha(args.campaign), 'Incomplete extension acceptance')
    need(core['accepted_tasks'] == 726 and core['lfine_valid_threshold_rows'] == 278784
         and core['gate_sha256'] == gate['release_gate']['sha256'], 'Incomplete core acceptance')
    need(sha(gate['task_manifest']['path']) == gate['task_manifest']['sha256'], 'Changed task roster')
    tasks = load(gate['task_manifest']['path'])
    need(len(tasks) == 2617, 'Incomplete task roster')
    need(args.aggregate == CAMP / 'evaluation/extensions' and args.selected == CAMP / 'evaluation/embedding_selected',
         'Unexpected prerequisite directories')
    inputs = {str(p): sha(p) for p in [args.campaign, extension_path, core_path,
        Path(gate['task_manifest']['path']), Path(gate['release_gate']['path'])]}
    aggregate = checked(args.aggregate, inputs, 'completed')
    need(aggregate['tasks'] == 2617 and aggregate['campaign_gate_sha256'] == sha(args.campaign)
         and aggregate['completion_sha256'] == sha(extension_path), 'Different extension aggregation')
    evaluation = checked(CAMP / 'evaluation', inputs, 'completed')
    need(evaluation['n_units_evaluated'] == 726 and evaluation['n_valid_threshold_rows'] == 278784
         and evaluation['frozen_semantics'] == aggregate['frozen_semantics'], 'Different/incomplete Lfine core')
    selected = checked(args.selected, inputs, 'passed_independent_selection')
    need(selected['n_representations'] == 1694 and selected['n_selection_rows'] == 420
         and selected['n_paired_contrasts'] == 84, 'A1 selection incomplete')
    for path in [args.campaign, args.aggregate / 'manifest.json', CAMP / 'evaluation/manifest.json']:
        need(selected['inputs'].get(str(path)) == sha(path), 'A1 selection belongs to different inputs')
    selected_nb = load(args.selected / 'notebook_receipt.json')
    notebook = ROOT / 'notebooks/dgscrna_results.ipynb'
    need(selected_nb['status'] == 'updated' and sha(notebook) == selected_nb['notebook_candidate_sha256'],
         'Notebook differs from the completed A1 update')
    previous_candidate = args.selected / 'dgscrna_results_candidate.ipynb'
    need(sha(previous_candidate) == selected_nb['notebook_candidate_sha256'], 'A1 immutable notebook differs')
    inputs[str(previous_candidate)] = sha(previous_candidate)
    metadata = load(notebook)['metadata']
    need(metadata['packaged_core_rerun']['gate_sha256'] == gate['release_gate']['sha256']
         and metadata['packaged_core_rerun']['evaluation_manifest_sha256'] == sha(CAMP / 'evaluation/manifest.json')
         and metadata['packaged_A1_selection']['summary_sha256'] == sha(args.selected / 'manifest.json'),
         'Notebook provenance does not match completed science')
    lock = load(CODE / 'REPORT_LOCK.json')
    for path, digest in lock['files'].items():
        need(sha(path) == digest, 'Frozen report source changed: ' + path)
        inputs[path] = digest
    import importlib.metadata
    need(str(Path(sys.executable).absolute()) == lock['runtime']['python'], 'Different report interpreter')
    for package, version in lock['runtime']['versions'].items():
        need(importlib.metadata.version(package) == version, 'Different report dependency: ' + package)
    source = load(CODE / 'SOURCE_MANIFEST.json')
    for path, digest in source['source_hashes'].items():
        need(sha(path) == digest, 'Old scientific source changed: ' + path)
        inputs[path] = digest
    for name, digest in source['derived_files'].items():
        need(sha(CODE / name) == digest, 'Derived source changed: ' + name)
    inputs.update({str(p): sha(p) for p in [CODE / 'REPORT_LOCK.json', CODE / 'SOURCE_MANIFEST.json',
        args.selected / 'notebook_receipt.json']})
    need(not (CAMP / 'evaluation/final_report').exists(), 'Preserve an existing final report attempt')
    return gate, tasks, aggregate, inputs


def old_inputs(inputs):
    a2 = OLD / 'no_clustering_lfine_v1'
    comp = OLD / 'comparison_lfine_v1'
    ap = load(a2 / 'full_validation/validation.json')
    cp = load(comp / 'selection_validation/validation.json')
    need(ap['status'] == cp['status'] == 'passed', 'Original independent inference validation missing')
    need(ap['summary_manifest_sha256'] == sha(a2 / 'summary/manifest.json')
         and cp['summary_manifest_sha256'] == sha(comp / 'summary/manifest.json'), 'Original inference proof mismatch')
    checked(a2 / 'summary', inputs, 'completed')
    checked(comp / 'summary', inputs, 'completed_pending_independent_validation')
    need(cp['n_fold_choices'] == 80 and cp['n_patient_endpoint_rows'] == 912, 'Original C validation roster differs')
    for proof, path in [(ap, a2 / 'full_validation/validation.json'), (cp, comp / 'selection_validation/validation.json')]:
        inputs[str(path)] = sha(path)
        for name, digest in proof.get('files', proof.get('outputs', {})).items():
            source = path.parent / name
            need(sha(source) == digest, 'Original independent verification artifact changed')
            inputs[str(source)] = digest
    return a2, comp


def read_table(path, round_trip=False):
    import pandas as pd
    kwargs = dict(dtype={'sample': str, 'patient': str, 'cutoff': str})
    if round_trip:
        kwargs['float_precision'] = 'round_trip'
    return pd.read_csv(path, **kwargs)


def core_rows(tasks, inputs, semantics, round_trip):
    import pandas as pd
    frames = []
    roster = sorted({(t['sample'], t['budget']) for t in tasks if t['family'] == 'no_cluster'})
    need(len(roster) == 242, 'Core anchor roster differs')
    for sample, budget in roster:
        if round_trip and budget != 'hvg2000':
            continue
        directory = CAMP / 'evaluation/units' / sample / budget
        receipt = checked(directory, inputs, 'completed')
        need(receipt['sample'] == sample and receipt['budget'] == budget and receipt['n_valid'] == 384
             and receipt['frozen_semantics'] == semantics, 'Core unit endpoint mismatch')
        frame = read_table(directory / 'metrics.csv.gz', round_trip)
        need(len(frame) == 384 and frame.terminal_valid.all(), 'Core unit incomplete')
        if round_trip:
            frame = frame[frame.stage.eq('terminal090') & frame.route.eq('UMAP2_HDBSCAN_R')]
        else:
            frame = frame[frame.library.eq('CM2_glioma_other') & frame.cutoff.eq('mean')]
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def a2_rows(tasks, aggregate, inputs):
    import pandas as pd
    frames = []
    for task in tasks:
        if task['family'] != 'no_cluster':
            continue
        directory = Path(task['output']) / 'evaluation'
        for name in ['manifest.json', 'metrics.csv.gz']:
            path = directory / name
            need(aggregate['input_hashes'].get(str(path)) == sha(path), 'A2 unit differs from accepted aggregate')
        receipt = checked(directory, inputs, 'completed')
        need(receipt['family'] == 'no_cluster_seed_control' and receipt['n_valid'] == 10
             and receipt['frozen_semantics'] == aggregate['frozen_semantics'], 'Incomplete A2 unit')
        frame = read_table(directory / 'metrics.csv.gz')
        need(len(frame) == 10 and frame.terminal_valid.all(), 'Invalid A2 endpoint')
        need(frame.cutoff.str.startswith('cellwise_lambda_').all(), 'A2 lambda encoding differs')
        frame['lambda_value'] = frame.cutoff.str.removeprefix('cellwise_lambda_').astype(float)
        need(set(frame.lambda_value) == {0, .5, 1, 1.5, 2}, 'A2 lambda roster changed')
        frames.append(frame)
    need(len(frames) == 242, 'Incomplete A2 unit roster')
    return pd.concat(frames, ignore_index=True)


def b_rows(data, tasks):
    frame = data[data.family.isin(['neighbor_control', 'learning_control'])].copy()
    need(len(frame) == 408 and frame.terminal_valid.all(), 'B endpoint roster incomplete')
    taskmap = {t['index']: t for t in tasks}
    for index, row in frame.iterrows():
        task = taskmap[int(row.campaign_task)]
        config = task['configuration']
        need(task['sample'] == row['sample'] and task['budget'] == row.budget, 'B task identity differs')
        if row.family == 'neighbor_control':
            need(config['name'] == row.configuration and task['family'] == 'neighbors', 'B neighbor identity differs')
            values = dict(task='neighbors', control=row.configuration, snn_k=config['snn_k'],
                umap_neighbors=config['umap_neighbors'], model_seed=42, epochs=10)
        else:
            need(task['family'] == 'learning', 'B learning identity differs')
            seed, epoch = row.configuration.split('_epoch')
            need(seed == 'seed' + str(config['model_seed']), 'B learning seed differs')
            values = dict(task='learning', control=row.configuration, snn_k=20, umap_neighbors=30,
                model_seed=config['model_seed'], epochs=int(epoch))
        for key, value in values.items():
            frame.loc[index, key] = value
    for name in ['model_seed', 'epochs', 'snn_k', 'umap_neighbors']:
        frame[name] = frame[name].astype(int)
    keys = ['task', 'sample', 'budget', 'route', 'control', 'library', 'cutoff', 'stage']
    frame['row_key'] = frame[keys].astype(str).agg('|'.join, axis=1)
    need(frame.row_key.is_unique, 'Duplicate B terminal endpoint')
    return frame


def snippets(descriptions, a2_stats, choices, comparison, out):
    link = lambda name: '../' + str((out / name).relative_to(ROOT))
    lower = int(descriptions.default_relation.eq('lower_than_an_alternative').sum())
    b = ('## Neighborhood and training sensitivity\n\n'
         f'Fresh packaged runs: three size-selected samples, HVG2000/5000, two fixed marker contexts; terminal confidence 0.90. '
         f'Defaults fell below a tested alternative in {lower}/48 descriptive sweeps. These controls do not establish a cohort optimum.\n\n'
         '![GBM neighborhood sensitivity](attachment:neighbor_lfine.png)\n\n'
         f'132 unique conditions; 12 shared defaults appear twice. [All checkpoint/seed results]({link("B/checkpoint_lfine_all_thresholds.csv")}) '
         'retain epochs 5/10/20/30, seeds 0/1/42 and both thresholds. No truth-based early stopping; '
         'TKU4163/HVG5000 has one known seed class, so its training curves do not validate multiclass biology.')
    a = ['## Cellwise seed replacement', '',
        'Fresh packaged terminal outputs reproduce the frozen patient vectors and fold choices exactly. '
        'Prespecified λ=1, glioma/other markers, confidence 0.90; 97 samples / 55 patients. Positive differences favor cellwise seeding.', '',
        '| Original route | HVG2000: ΔF1 (95% CI), Holm p | HVG5000: ΔF1 (95% CI), Holm p |', '|---|---:|---:|']
    names = {'PCA30_HDBSCAN_R': 'PCA–HDBSCAN', 'PCA30_SNN': 'PCA–SNN', 'UMAP2_HDBSCAN_R': 'UMAP–HDBSCAN', 'UMAP2_SNN': 'UMAP–SNN'}
    fixed = a2_stats[a2_stats.cohort.eq('primary') & a2_stats.selection.eq('fixed_lambda1')]
    for route, label in names.items():
        values = []
        for budget in ['hvg2000', 'hvg5000']:
            part = fixed[fixed.reference_route.eq(route) & fixed.budget.eq(budget)]
            need(len(part) == 1, 'A2 display contrast missing')
            r = part.iloc[0]
            values.append(f'{r.mean_delta_cellwise_minus_original:+.4f} ({r.CI025:+.4f}, {r.CI975:+.4f}); {r.p_holm:.3g}')
        a.append('| ' + ' | '.join([label] + values) + ' |')
    selections = ', '.join(f'λ={value:g}: {count}/20 folds' for value, count in choices.lambda_value.value_counts().sort_index().items())
    a += ['', f'Training-patient selection is separate ({selections}); valid no-training terminal states remain included. '
        f'[All 32 contrasts]({link("A2_retained_original_inference.csv")}) · [Fresh execution counts]({link("A2_fresh_execution_states.csv")}). '
        'Original CI/p are retained only after exact inference-input revalidation: 2,000 paired patient bootstrap draws, '
        'Pratt Wilcoxon, Holm16 per cohort; conditional on saved predictions/configurations.']
    c = ['## Method comparison · original HVG2000 partition', '',
        '97 samples / 55 patients; five patient folds. Cluster tools share HVG2000 → PCA30 → UMAP2 → HDBSCAN. '
        'DG-scRNA uses fresh packaged terminal outputs; the five comparator methods retain their archived fixed predictions. '
        'Frozen DG marker/cutoff choices and held-out patient vectors reproduce exactly.', '',
        '| Method | Lfine F1 | Coverage | Off-vocab | ΔF1 [95% CI] | Holm p |', '|---|---:|---:|---:|---:|---:|']
    for row in comparison:
        delta = '—' if row['method'] == 'DG-scRNA' else f"{row['mean_delta']:+.3f} [{row['CI95_low']:+.3f}, {row['CI95_high']:+.3f}]"
        p = '—' if row['method'] == 'DG-scRNA' else f"{row['p_Holm']:.3g}"
        c.append(f"| {row['method']} | {row['lfine_macroF1']:.3f} | {row['coverage']:.1%} | {row['offvocab_rate']:.1%} | {delta} | {p} |")
    c += ['', 'ΔF1 = competitor − DG-scRNA. Coverage includes non-abstaining unmappable calls; Off-vocab reports those calls separately. '
        'Scores measure Lfine compatible-target-set performance. The marker methods share 13 libraries; SingleR uses labelled training patients; '
        'scDeepSort uses a pretrained Brain atlas without malignant classes.', '',
        'Original CI/p are retained after exact DG patient-vector revalidation: 10,000 paired bootstrap draws, default Wilcoxon, Holm5; '
        'selection is not repeated within resampling. These comparisons do not establish overall superiority. '
        f'[Current DG and archived comparator table]({link("C_current_DG_archived_comparators.csv")}). '
        'The linked earlier equal24 analysis remains a separate archived sensitivity analysis.']
    robustness = ('# 6 · Robustness\n\n'
        'Complete fresh package campaign: 726 core units and 2,617 reviewer-control units passed terminal-output parity and Lfine evaluation. '
        'The primary cohort remains 97 samples / 55 patients; all 121 samples and the 24-sample sensitivity subset are retained. '
        'The latter is not an independent patient cohort. §9 reports the full feature/route grid and completed reviewer comparisons; '
        'non-DG comparator predictions retain their archived provenance.')
    return {'# 6 · Robustness': robustness, '## Neighborhood and training sensitivity': b,
        '## Cellwise seed replacement': '\n'.join(a), '## Method comparison · original HVG2000 partition': '\n'.join(c)}


def update_notebook(out, texts, update, expected_before=None):
    import nbformat
    notebook = ROOT / 'notebooks/dgscrna_results.ipynb'
    lockpath = OLD / 'notebook/update.lock'
    with lockpath.open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        before = sha(notebook)
        if expected_before is not None:
            need(before == expected_before, 'Notebook changed during report revalidation; review the newer edit')
        report = nbformat.read(notebook, as_version=4)
        need(len(report.cells) == 20 and report.metadata.get('compact_reference_report'), 'Unexpected notebook structure')
        retained_metadata = {key: json.dumps(report.metadata[key], sort_keys=True)
            for key in ['packaged_core_rerun', 'packaged_A1_selection']}
        changes = []
        untouched = {i: json.dumps(cell, sort_keys=True) for i, cell in enumerate(report.cells)}
        for prefix, text in texts.items():
            matches = [i for i, cell in enumerate(report.cells) if cell.source.startswith(prefix)]
            need(len(matches) == 1, 'Existing report section missing/duplicated: ' + prefix)
            index = matches[0]
            report.cells[index].source = text
            if prefix == '## Neighborhood and training sensitivity':
                report.cells[index]['attachments'] = {'neighbor_lfine.png': {
                    'image/png': base64.b64encode((out / 'B/neighbor_lfine.png').read_bytes()).decode()}}
            changes.append(index)
        for index, value in untouched.items():
            if index not in changes:
                need(json.dumps(report.cells[index], sort_keys=True) == value, 'Unrelated notebook cell changed')
        for key, value in retained_metadata.items():
            need(json.dumps(report.metadata[key], sort_keys=True) == value, 'Core/A1 metadata changed')
        report.metadata['packaged_extension_report'] = dict(summary_sha256=sha(out / 'manifest.json'),
            B_fresh_terminal_rows=408, A2_exact_patient_vectors=True, C_exact_DG_patient_vectors=True,
            non_DG_comparators_refitted=False, original_inference_retained=True)
        nbformat.validate(report)
        candidate = out / 'dgscrna_results_candidate.ipynb'
        nbformat.write(report, candidate)
        need(sha(notebook) == before, 'Concurrent notebook modification')
        if update:
            shutil.copy2(notebook, out / ('notebook_before_' + before + '.ipynb'))
            temporary = notebook.with_name('.dgscrna_results.final_report.tmp')
            shutil.copyfile(candidate, temporary)
            temporary.replace(notebook)
        write(out / 'notebook_receipt.json', dict(status='updated' if update else 'candidate',
            notebook_before_sha256=before, notebook_candidate_sha256=sha(candidate), replaced_cells=changes,
            remaining_cells_preserved=True, core_A1_metadata_preserved=True, PTC_changed=False,
            HTML_created=False, OneDrive_accessed=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', type=Path, required=True)
    parser.add_argument('--aggregate', type=Path, required=True)
    parser.add_argument('--selected', type=Path, required=True)
    parser.add_argument('--update-notebook', action='store_true')
    args = parser.parse_args()
    for name in ['campaign', 'aggregate', 'selected']:
        setattr(args, name, getattr(args, name).resolve())
    gate, tasks, aggregate, inputs = prerequisites(args)
    a2, comp = old_inputs(inputs)
    import pandas as pd
    from render_B import render
    from revalidate import revalidate_a2, revalidate_c
    out = CAMP / 'evaluation/final_report'
    out.mkdir(parents=True)
    data = read_table(args.aggregate / 'metrics_all_extensions.csv.gz')
    b = b_rows(data, tasks)
    b.to_csv(out / 'B_fresh_terminal_metrics.csv.gz', index=False)
    descriptions, checkpoints = render(b, out / 'B')
    readold = lambda directory, name: read_table(directory / 'summary' / name, round_trip=True)
    a2_candidates = readold(a2, 'all_candidate_metrics.csv')
    a2_candidates = a2_candidates[a2_candidates.stage.isin(STAGES)].copy()
    a2_anchors = readold(a2, 'all_original_route_metrics.csv')
    a2_anchors = a2_anchors[a2_anchors.stage.isin(STAGES)].copy()
    a2_stats = readold(a2, 'patient_paired_summary.csv')
    a2_proof, choices = revalidate_a2(a2_rows(tasks, aggregate, inputs),
        core_rows(tasks, inputs, aggregate['frozen_semantics'], False), a2_candidates, a2_anchors,
        readold(a2, 'training_patient_lambda_choices.csv'), readold(a2, 'patient_paired_results.csv'),
        a2_stats, read_table(OLD / 'no_clustering/patient_folds.csv'), out)
    c_proof, comparison = revalidate_c(core_rows(tasks, inputs, aggregate['frozen_semantics'], True),
        readold(comp, 'all_eligible_candidate_metrics.csv.gz'), readold(comp, 'training_patient_choices.csv'),
        readold(comp, 'patient_heldout_results.csv'), readold(comp, 'patient_heldout_summary.csv'),
        readold(comp, 'paired_patient_comparisons.csv'),
        read_table(ROOT / 'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/protocol/patient_folds.csv', True), out)
    texts = snippets(descriptions, a2_stats, choices, comparison, out)
    write(out / 'section_texts.json', texts)
    proof = dict(B=dict(fresh_rows=len(b), unique_neighbor_conditions=132, displayed_neighbor_points=144,
        descriptive_sweeps=len(descriptions), checkpoint_rows=len(checkpoints)), A2=a2_proof, C=c_proof)
    write(out / 'numerical_revalidation.json', proof)
    for path, digest in inputs.items():
        need(sha(path) == digest, 'Input changed during report generation: ' + path)
    sources = {str(path): sha(path) for path in CODE.glob('*') if path.is_file()}
    value = dict(status='passed_existing_sections_revalidated', campaign_gate_sha256=sha(args.campaign),
        core_units=726, extension_units=2617, numerical_revalidation=proof,
        source_hashes=sources, input_hashes=inputs,
        inference_policy='Original A2/C CI and p values retained only after exact float64 patient vectors and fold choices match; non-DG predictions remain archived',
        parser_policy='A2 original per-unit default CSV parse; C original per-unit round_trip parse and ordered reductions; generic aggregate required and bound',
        outputs={str(p.relative_to(out)): sha(p) for p in out.rglob('*') if p.is_file()},
        job=os.environ['SLURM_JOB_ID'], step=os.environ.get('SLURM_STEP_ID'), created_at=datetime.now(timezone.utc).isoformat())
    write(out / 'manifest.json', value)
    update_notebook(out, texts, args.update_notebook,
        expected_before=inputs[str(args.selected / 'dgscrna_results_candidate.ipynb')])
    (out / 'COMPLETE').write_text(sha(out / 'manifest.json') + '\n')
    print(json.dumps(dict(status=value['status'], numerical_revalidation=proof)))


if __name__ == '__main__':
    main()
