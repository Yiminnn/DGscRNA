#!/usr/bin/env python3
"""SLURM-only Lfine evaluation of freshly packaged GBM terminal endpoints.

Uses frozen v5 compatible-target-set semantics, not strict fine-type accuracy.
`unit` evaluates one core task; `run` evaluates a completed selected-arm pilot;
`aggregate` includes every prespecified core condition, retaining unavailable NA.
This adapter never loads expression matrices, fits models, or changes predictions.
"""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import itertools
import json
import os
import tempfile

ROOT = Path(__file__).resolve().parents[2]
NATIVE = ROOT / 'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917'
CAMPAIGN = ROOT / 'results/hvg_ptc_20260916_v1/package_reference_rerun_20260920'
TASKS = ROOT / 'handoff/package_release_20260920/rerun_inventory.core_tasks.json'
SOURCE_LOCK = Path(__file__).with_name('evaluate_core_lfine_sources.json')
ROUTES = ['PCA30_SNN', 'PCA30_HDBSCAN_R', 'UMAP2_SNN', 'UMAP2_HDBSCAN_R']
CUTS = ['none', 'mean', '0.5']
STAGES = [('terminal070', 'final070', .70), ('terminal090', 'final090', .90)]
KEYS = ['sample', 'budget', 'route', 'library', 'cutoff', 'stage']
GROUPS = ['budget', 'route', 'library', 'cutoff', 'stage']
FIXED = 'CM2_glioma_other'


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write_json(path, value):
    path = Path(path)
    with tempfile.NamedTemporaryFile('w', dir=path.parent, delete=False) as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        temporary = Path(stream.name)
    temporary.replace(path)


def checked_json(directory, manifest, flag):
    directory = Path(directory)
    path = directory / manifest
    digest = sha(path)
    if (directory / flag).read_text().strip() != digest:
        raise ValueError(f'Completion checksum mismatch: {path}')
    return json.loads(path.read_text()), digest


def initialize():
    """Import only the old side-effect-free provider; pin every semantic input."""
    global np, pd, compact, v5, helpers, mappings, libraries, cohort, reference_overlap
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Scientific evaluation requires SLURM')
    import numpy as np
    import pandas as pd
    lock = json.loads(SOURCE_LOCK.read_text())
    for relative, digest in lock['sources'].items():
        if sha(ROOT / relative) != digest:
            raise ValueError(f'Frozen evaluation source changed: {relative}')
    source = ROOT / 'handoff/lfine_compact_20260920/evaluate_lfine.py'
    spec = importlib.util.spec_from_file_location('frozen_compact_lfine', source)
    compact = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(compact)
    compact.initialize()
    v5, helpers, mappings, libraries = compact.v5, compact.helpers, compact.mappings, compact.libraries
    marker_manifest = json.loads((NATIVE / 'markers' / 'manifest.json').read_text())
    reference_overlap = {source['library']: source.get('reference_label_construction_overlap')
        for source in marker_manifest['sources']}
    assert set(reference_overlap) == set(libraries)
    cohort = pd.read_csv(NATIVE / 'protocol/cohort.csv', dtype={'sample': str, 'patient': str})
    assert len(cohort) == 121 and int(cohort.primary.sum()) == 97
    assert cohort['sample'].is_unique and len(libraries) == 16
    return lock


def get_task(index):
    task = json.loads(TASKS.read_text())[index]
    assert task['index'] == index
    return task


def expected_rows(task, routes=None, selected_libraries=None, cutoffs=None):
    member = cohort[cohort['sample'].eq(task['sample'])]
    assert len(member) == 1
    assert member.iloc[0].patient == task['patient']
    assert bool(member.iloc[0].primary) == task['primary']
    selected_libraries = libraries if selected_libraries is None else selected_libraries
    result = []
    for route, library, cutoff, (stage, _, threshold) in itertools.product(
            routes or task['routes'], selected_libraries, cutoffs or task['cutoffs'], STAGES):
        result.append(dict(sample=task['sample'], patient=task['patient'], primary=task['primary'],
            task_index=task['index'], budget=task['budget'], route=route, library=library,
            cutoff=cutoff, stage=stage, terminal_threshold=threshold, n_cells=task['n_cells'],
            family='native_R_budget', status='not_evaluated', terminal_valid=False,
            dl_status='not_recorded', training_executed=False, n_unique_models_trained=0,
            identical_result_reused=False, unavailable_reason='missing_evaluation_receipt',
            core_hvg24=(library == FIXED and cutoff == 'mean'),
            fullgene_marker16=(task['budget'] == 'all' and route == 'UMAP2_HDBSCAN_R' and cutoff == 'mean'),
            reference_label_construction_overlap=reference_overlap[library],
            **{name: np.nan for name in compact.SCORE_FIELDS}))
    return result


def truth_for(task):
    truth_path = NATIVE / 'evaluation_inputs' / task['sample'] / 'truth.csv.gz'
    manifest = json.loads((NATIVE / 'inputs' / task['sample'] / 'input_manifest.json').read_text())
    digest = sha(truth_path)
    assert digest == manifest['evaluation_files']['truth.csv.gz']
    truth = pd.read_csv(truth_path, dtype=str, keep_default_na=False).rename(
        columns={'cell_id': 'CellID', 'lfine_original': 'Lfine'})
    assert truth.CellID.is_unique and truth.Patient.eq(task['patient']).all()
    assert len(truth) == task['n_cells']
    scope = v5.label_scope(truth, helpers)
    assert np.array_equal(scope[0], truth.Lfine.to_numpy())
    return truth, scope, digest


def score_predictions(prediction, library, truth, scope):
    mapped = np.asarray([mappings[library].get(value,
        'Unknown' if value in compact.ABSTAIN else 'UNMAPPABLE') for value in prediction], dtype=object)
    values = v5.metrics(mapped, truth, scope, helpers, library)
    # v5 also computes a legacy coarse malignant metric internally. Never emit it.
    result = {name: values[name] for name in compact.SCORE_FIELDS}
    assert result['n_reference_lfine_disagreements'] == 0
    return result


def evaluate(task, run_root, selected=False, compare_reference=False):
    run_root = Path(run_root).resolve()
    kwargs = {}
    if selected:
        cfg = json.loads((run_root / 'run_config.json').read_text())['config']
        assert cfg['sample'] == task['sample']
        kwargs = dict(routes=ROUTES if cfg['route'] == 'all' else [cfg['route']],
            selected_libraries=libraries if cfg['library'] == 'all' else [cfg['library']],
            cutoffs=CUTS if cfg['cutoff'] == 'all' else [cfg['cutoff']])
    rows = expected_rows(task, **kwargs)
    truth, scope, truth_hash = truth_for(task)
    for row in rows:
        row.update(run_root=str(run_root), truth_sha256=truth_hash, lfine_n_classes=len(scope[1]))
    if not (run_root / 'COMPLETE').exists():
        for row in rows:
            row.update(status='unavailable', unavailable_reason='run_not_complete')
        return rows, []
    run, run_hash = checked_json(run_root, 'run_manifest.json', 'COMPLETE')
    assert run['status'] == 'completed' and run['reference_labels_used_for_fit'] is False
    config = json.loads((run_root / 'run_config.json').read_text())['config']
    assert config['sample'] == task['sample']
    assert ('all' if config['features'] == 'all' else 'hvg' + config['features']) == task['budget']
    actual = {(c['route'], c['library'], str(c['cutoff'])): c for c in run['conditions']}
    requested = {(r['route'], r['library'], r['cutoff']) for r in rows}
    assert len(actual) == len(run['conditions']) == run['terminal_condition_count']
    assert set(actual) == requested, 'Completed run arm set differs from prespecified evaluation arm set'
    cache, old_cache, parity = {}, {}, []
    for row in rows:
        route, library, cutoff = row['route'], row['library'], row['cutoff']
        source = run_root / 'GBM' / task['sample'] / task['budget'] / route
        key = (route, library, cutoff)
        row.update(run_manifest_sha256=run_hash, package_version=run['package_version'])
        try:
            if route not in cache:
                score, sh = checked_json(source, 'score_manifest.json', 'SCORE_COMPLETE')
                cache[route] = (score, sh)
            score, sh = cache[route]
            ids = [aid for aid, arm in score['arms'].items()
                   if arm['library'] == library and str(arm['cutoff']) == cutoff]
            assert len(ids) == 1
            aid = ids[0]
            dest = source / 'terminal' / aid
            row.update(arm_id=aid, score_manifest_sha256=sh)
            if key not in cache:
                terminal, th = checked_json(dest, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
                # Preserve declared execution status even when later artifact checks fail.
                row.update(terminal_manifest_sha256=th, dl_status=terminal.get('dl_status', 'not_recorded'),
                    training_executed=terminal.get('training_executed') is True,
                    terminal_declared_valid=terminal.get('terminal_valid') is True)
                assert terminal['score_manifest_sha256'] == sh
                assert terminal['arm']['library'] == library and str(terminal['arm']['cutoff']) == cutoff
                assert terminal['reference_labels_used_for_fit'] is False
                if terminal['status'] != 'completed' or not terminal['terminal_valid']:
                    raise ValueError('invalid_terminal_state')
                trained_states = {'trained', 'trained_single_known_class'}
                terminal_states = trained_states | {'no_op_all_initially_known',
                    'no_known_labels_archived_Undecided_terminal', 'structural_insufficient_known_split'}
                if terminal['dl_status'] not in terminal_states:
                    raise ValueError(f"invalid_terminal_execution:{terminal['dl_status']}")
                assert terminal['training_executed'] == (terminal['dl_status'] in trained_states)
                if terminal['dl_status'] == 'no_op_all_initially_known':
                    assert terminal['n_pool'] == 0
                elif terminal['dl_status'] == 'no_known_labels_archived_Undecided_terminal':
                    assert terminal['n_known'] == 0
                elif terminal['dl_status'] == 'structural_insufficient_known_split':
                    assert terminal['n_known'] == 1
                pp = dest / 'predictions.csv.gz'
                ph = sha(pp)
                assert ph == terminal['predictions_sha256'] == actual[key]['predictions_sha256']
                assert (run_root / actual[key]['predictions']).resolve() == pp.resolve()
                pred = pd.read_csv(pp, dtype=str, keep_default_na=False,
                    usecols=['cell_id', 'initial', 'final070', 'final090'])
                assert np.array_equal(pred.cell_id, truth.CellID)
                known = pred.initial.ne('Undecided')
                for column in ['final070', 'final090']:
                    assert np.array_equal(pred.loc[known, 'initial'], pred.loc[known, column])
                cache[key] = (terminal, th, pred, ph, pp)
            terminal, th, pred, ph, pp = cache[key]
            stage_column = 'final070' if row['stage'] == 'terminal070' else 'final090'
            values = score_predictions(pred[stage_column], library, truth, scope)
            row.update(**values, status='completed', terminal_valid=True, unavailable_reason='',
                terminal_manifest_sha256=th, terminal_declared_valid=True,
                predictions_sha256=ph, predictions_path=str(pp),
                dl_status=terminal['dl_status'], training_executed=terminal['training_executed'],
                identical_result_reused=terminal['identical_result_reused'],
                n_unique_models_trained=int(terminal['training_executed'] and not terminal['identical_result_reused']),
                DL_features=terminal['DL_features'], n_training_classes=terminal['n_training_classes'],
                n_pool=terminal['n_pool'], source_native_labels_preserved=True)
            if compare_reference:
                old = Path(task['expected_reference']) / route
                if route not in old_cache:
                    old_cache[route] = checked_json(old, 'score_manifest.json', 'SCORE_COMPLETE')[0]
                old_ids = [k for k, arm in old_cache[route]['arms'].items()
                    if arm['library'] == library and str(arm['cutoff']) == cutoff]
                assert len(old_ids) == 1
                old_terminal = old / 'terminal' / old_ids[0]
                old_manifest, _ = checked_json(old_terminal, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
                old_pp = old_terminal / 'predictions.csv.gz'
                assert sha(old_pp) == old_manifest['predictions_sha256']
                old_pred = pd.read_csv(old_pp, dtype=str, keep_default_na=False,
                    usecols=['cell_id', stage_column])
                assert np.array_equal(old_pred.cell_id, truth.CellID)
                assert np.array_equal(old_pred[stage_column], pred[stage_column])
                old_values = score_predictions(old_pred[stage_column], library, truth, scope)
                np.testing.assert_allclose([old_values[k] for k in compact.SCORE_FIELDS],
                    [values[k] for k in compact.SCORE_FIELDS], rtol=0, atol=0, equal_nan=True)
                parity.append({k: row[k] for k in KEYS} | dict(
                    reference_terminal_labels_exact=True, reference_lfine_metrics_exact=True,
                    n_metric_checks=len(compact.SCORE_FIELDS)))
        except (OSError, ValueError, AssertionError, KeyError) as exc:
            # Never silently substitute marker calls for missing or invalid terminal calls.
            row.update(status='invalid', terminal_valid=False, unavailable_reason=f'{type(exc).__name__}:{exc}')
            row.update({name: np.nan for name in compact.SCORE_FIELDS})
    return rows, parity


def save_evaluation(out, rows, parity, lock, task, compare_reference=False):
    out.mkdir(parents=True, exist_ok=True)
    (out / 'COMPLETE').unlink(missing_ok=True)
    frame = pd.DataFrame(rows)
    assert not frame.duplicated(KEYS).any()
    frame.to_csv(out / 'metrics.csv.gz', index=False)
    if compare_reference:
        pd.DataFrame(parity).to_csv(out / 'reference_parity.csv', index=False)
    outputs = {'metrics.csv.gz': sha(out / 'metrics.csv.gz')}
    if compare_reference:
        outputs['reference_parity.csv'] = sha(out / 'reference_parity.csv')
    manifest = dict(status='completed' if frame.terminal_valid.all() else 'completed_with_unavailable_conditions',
        task_index=task['index'], sample=task['sample'], budget=task['budget'],
        n_conditions=len(frame), n_valid=int(frame.terminal_valid.sum()),
        status_counts=frame.status.value_counts().to_dict(),
        dl_status_counts=frame.dl_status.value_counts().to_dict(),
        reference_parity_conditions=len(parity), reference_parity_requested=compare_reference,
        frozen_semantics=lock, script_sha256=sha(__file__), outputs=outputs,
        no_fitting=True, no_expression_matrix_loading=True, no_truth_based_model_or_marker_selection=True,
        display_policy='Lfine only; semantic-parent source fields are ontology provenance, not reported coarse scores',
        timestamp=datetime.now(timezone.utc).isoformat(), job=os.environ['SLURM_JOB_ID'])
    write_json(out / 'manifest.json', manifest)
    (out / 'COMPLETE').write_text(sha(out / 'manifest.json') + '\n')
    return manifest


def summarize(data, keys):
    # Same sample and patient averaging as compact.summarize, with status breakdowns.
    result = compact.summarize(data, keys)
    return result


def aggregate(out, lock):
    tasks = json.loads(TASKS.read_text())
    frames, audits = [], []
    for task in tasks:
        directory = out / 'units' / task['sample'] / task['budget']
        if (directory / 'COMPLETE').exists():
            manifest, _ = checked_json(directory, 'manifest.json', 'COMPLETE')
            assert manifest['frozen_semantics'] == lock
            assert manifest['script_sha256'] == sha(__file__), 'Re-evaluate units after adapter changes'
            assert manifest['sample'] == task['sample'] and manifest['budget'] == task['budget']
            assert sha(directory / 'metrics.csv.gz') == manifest['outputs']['metrics.csv.gz']
            frame = pd.read_csv(directory / 'metrics.csv.gz', dtype={'cutoff': str}, keep_default_na=True)
            assert len(frame) == 384 and not frame.duplicated(KEYS).any()
            expected = pd.DataFrame(expected_rows(task))
            assert set(map(tuple, frame[KEYS].values)) == set(map(tuple, expected[KEYS].values))
            audits.append(dict(task_index=task['index'], sample=task['sample'], budget=task['budget'],
                status=manifest['status'], n_valid=manifest['n_valid']))
        else:
            frame = pd.DataFrame(expected_rows(task))
            audits.append(dict(task_index=task['index'], sample=task['sample'], budget=task['budget'],
                status='not_evaluated', n_valid=0))
        frames.append(frame)
    data = pd.concat(frames, ignore_index=True)
    assert len(data) == 278784 and not data.duplicated(KEYS).any()
    assert data.terminal_valid.dtype == bool and data.primary.dtype == bool
    assert len(data['sample'].unique()) == 121
    out.mkdir(parents=True, exist_ok=True)
    (out / 'COMPLETE').unlink(missing_ok=True)
    data.to_csv(out / 'metrics_all_conditions.csv.gz', index=False)
    summarize(data, GROUPS).to_csv(out / 'summary_all_conditions.csv', index=False)
    for name, mask, grouping in [
            ('hvg24', data.core_hvg24, ['budget', 'route', 'stage']),
            ('fullgene_markers', data.fullgene_marker16, ['library', 'stage'])]:
        subset = data[mask].copy()
        subset.to_csv(out / f'metrics_{name}.csv.gz', index=False)
        summarize(subset, grouping).to_csv(out / f'summary_{name}.csv', index=False)
    pd.DataFrame(audits).to_csv(out / 'task_status.csv', index=False)
    status = data.groupby(GROUPS + ['status', 'dl_status'], dropna=False).size().rename('n_samples').reset_index()
    status.to_csv(out / 'status_counts.csv', index=False)
    compact.OUT = out
    mapping_audit = compact.audit_mapping()
    historical_checks = compact.verify_historical()
    manifest = dict(status='completed' if data.terminal_valid.all() else 'partial',
        n_expected_units=726, n_expected_terminal_conditions=139392, n_expected_threshold_rows=len(data),
        n_valid_threshold_rows=int(data.terminal_valid.sum()), n_units_evaluated=sum(a['status'] != 'not_evaluated' for a in audits),
        n_samples=121, n_primary_samples=97, n_sensitivity_samples=24,
        endpoint='v5 compatible-target-set Lfine macro-F1; observed classes support>=20 except Other/nan; all cells remain in TP/FP/FN',
        cohort_aggregation='Equal patient weights after within-patient available-sample means; fixed historical 97/24 split; sample means also supplied',
        incomplete_cohort_policy='Expected, valid and unavailable counts reported per condition; partial comparisons do not establish a cohort winner',
        display_policy='No author coarse-label performance metric is emitted',
        invalid_terminal_policy='Missing/invalid terminal results remain NA, never marker-only substitutes',
        thresholds='Both 0.70 and 0.90 derive from one terminal model; they are not separate fits',
        marker_dependency='Frozen source provenance flags BrainAtlas112, CARE_TME and UNION_all for author-reference overlap; unspecified overlap remains unknown, not independent validation',
        historical_v5_checks=historical_checks, mapping_audit=mapping_audit,
        frozen_semantics=lock, script_sha256=sha(__file__),
        outputs={p.name: sha(p) for p in out.iterdir() if p.is_file() and p.name not in {'manifest.json', 'COMPLETE'}},
        timestamp=datetime.now(timezone.utc).isoformat(), job=os.environ['SLURM_JOB_ID'])
    write_json(out / 'manifest.json', manifest)
    (out / 'COMPLETE').write_text(sha(out / 'manifest.json') + '\n')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    unit = sub.add_parser('unit')
    unit.add_argument('--task-index', type=int, required=True)
    unit.add_argument('--campaign-root', type=Path, default=CAMPAIGN)
    unit.add_argument('--out', type=Path)
    unit.add_argument('--compare-reference', action='store_true')
    run = sub.add_parser('run')
    run.add_argument('--run-root', type=Path, required=True)
    run.add_argument('--out', type=Path, required=True)
    run.add_argument('--compare-reference', action='store_true')
    total = sub.add_parser('aggregate')
    total.add_argument('--out', type=Path, default=CAMPAIGN / 'evaluation')
    args = parser.parse_args()
    lock = initialize()
    if args.command == 'aggregate':
        manifest = aggregate(args.out.resolve(), lock)
    else:
        if args.command == 'unit':
            task = get_task(args.task_index)
            run_root = args.campaign_root / 'core' / task['sample'] / task['budget']
            out = args.out or args.campaign_root / 'evaluation' / 'units' / task['sample'] / task['budget']
        else:
            run_root, out = args.run_root, args.out
            cfg = json.loads((run_root / 'run_config.json').read_text())['config']
            budget = 'all' if cfg['features'] == 'all' else 'hvg' + cfg['features']
            matches = [t for t in json.loads(TASKS.read_text()) if t['sample'] == cfg['sample'] and t['budget'] == budget]
            assert len(matches) == 1
            task = matches[0]
        rows, parity = evaluate(task, run_root, selected=args.command == 'run', compare_reference=args.compare_reference)
        manifest = save_evaluation(out.resolve(), rows, parity, lock, task, args.compare_reference)
    for relative, digest in lock['sources'].items():
        assert sha(ROOT / relative) == digest, f'Frozen evaluation input changed during scoring: {relative}'
    print(json.dumps({k: v for k, v in manifest.items() if k in {
        'status', 'n_conditions', 'n_valid', 'reference_parity_conditions', 'n_expected_threshold_rows',
        'n_valid_threshold_rows', 'n_units_evaluated', 'historical_v5_checks'}}, indent=2), flush=True)
    if manifest['status'] == 'completed_with_unavailable_conditions':
        raise SystemExit(2)


if __name__ == '__main__':
    main()
