"""Lfine-only scoring for verified packaged GBM research controls.

The JSON specification contains sample, budget, family, configuration, and a
conditions list. Each condition names route/library/cutoff/terminal_directory.
An accepted parity receipt and source hashes must accompany the specification.
This evaluator never loads expression data, fits, selects, or relabels a model.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import importlib.util
import json
import os
import tempfile

HERE = Path(__file__).resolve().parent
CORE = HERE / 'evaluate_core_lfine.py'
CORE_SHA256 = '5352acee29aa5fb6e8f1e73549ea88e2c0801e037f7460335b4baadfcb172def'
LOCK_SHA256 = '8c1ccab90c51564e01e7547b3557933604c0cb263d4207db3d686867379a6712'
TRAINED = {'trained', 'trained_single_known_class'}
NOOP = {'no_op_all_initially_known', 'no_known_labels_archived_Undecided_terminal',
        'structural_insufficient_known_split'}


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def checked(directory, name, flag):
    directory = Path(directory)
    digest = sha(directory / name)
    require((directory / flag).read_text().strip() == digest,
            f'Invalid completion receipt: {directory}')
    return load(directory / name), digest


def evaluate(spec_path, output):
    require(os.environ.get('SLURM_JOB_ID'), 'Evaluation requires SLURM')
    evaluator_hash = sha(__file__)
    require(sha(CORE) == CORE_SHA256 and
            sha(HERE / 'evaluate_core_lfine_sources.json') == LOCK_SHA256,
            'Frozen core evaluation provider changed')
    spec_path, output = Path(spec_path).resolve(), Path(output).resolve()
    spec_hash = sha(spec_path)
    spec = load(spec_path)
    require(spec['family'] in {'geometry_only_fixed_DL2000', 'MLP_only',
            'representation_control', 'neighbor_control', 'learning_control',
            'no_cluster_seed_control', 'seven_space_control'}, 'Unspecified control family')
    proof_path = Path(spec['parity_receipt']['path']).resolve()
    require(sha(proof_path) == spec['parity_receipt']['sha256'], 'Changed parity receipt')
    proof = load(proof_path)
    require(proof['status'] in {'passed', 'passed_exact'}, 'Parity did not pass')
    require(proof['sample'] == spec['sample'] and proof['budget'] == spec['budget'],
            'Parity receipt belongs to a different sample/budget')
    for path, digest in spec['source_hashes'].items():
        require(sha(path) == digest, f'Changed control source: {path}')
    require(not output.exists(), 'Preserve existing evaluation; use a new output directory')
    conditions = spec['conditions']
    require(len(conditions) == spec['expected_conditions'] and len(conditions) > 0,
            'Incomplete expected condition roster')
    keys = [(c['route'], c['library'], str(c['cutoff'])) for c in conditions]
    require(len(set(keys)) == len(keys), 'Duplicate condition identities')
    parity_conditions = {(c['route'], c['library'], str(c['cutoff'])): c
                         for c in proof['conditions']}
    require(len(parity_conditions) == len(proof['conditions']), 'Duplicate parity condition identities')
    require(set(keys) <= set(parity_conditions) if spec.get('allow_proof_subset') is True
            else set(keys) == set(parity_conditions), 'Parity condition roster differs from evaluation')
    module_spec = importlib.util.spec_from_file_location('package_frozen_lfine', CORE)
    provider = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(provider)
    lock = provider.initialize()
    pd, np = provider.pd, provider.np
    tasks = load(provider.TASKS)
    candidates = [t for t in tasks if t['sample'] == spec['sample'] and t['budget'] == spec['budget']]
    require(len(candidates) == 1, 'Unknown sample/budget')
    task = candidates[0]
    truth, scope, truth_hash = provider.truth_for(task)
    rows, receipts = [], {}
    for condition in conditions:
        terminal_path = Path(condition['terminal_directory']).resolve()
        terminal, terminal_hash = checked(terminal_path, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
        parity_condition = parity_conditions[(condition['route'], condition['library'], str(condition['cutoff']))]
        require(parity_condition['status'] in {'passed', 'passed_exact'} and
                parity_condition['actual_manifest_sha256'] == terminal_hash,
                'Terminal result differs from the condition verified by parity')
        require(terminal['status'] == 'completed' and terminal['terminal_valid'] is True,
                'Invalid terminal result')
        require(terminal['reference_labels_used_for_fit'] is False, 'Reference labels used in fitting')
        require(terminal['dl_status'] in TRAINED | NOOP and
                terminal['training_executed'] == (terminal['dl_status'] in TRAINED),
                'Inconsistent DL execution state')
        require(terminal['arm']['library'] == condition['library'] and
                str(terminal['arm']['cutoff']) == str(condition['cutoff']),
                'Terminal marker condition changed')
        prediction_path = terminal_path / 'predictions.csv.gz'
        require(sha(prediction_path) == terminal['predictions_sha256'], 'Prediction checksum mismatch')
        calls = pd.read_csv(prediction_path, dtype=str, keep_default_na=False,
                           usecols=['cell_id', 'initial', 'final070', 'final090'])
        require(np.array_equal(calls.cell_id, truth.CellID), 'Cell identity/order mismatch')
        known = calls.initial.ne('Undecided')
        require(int((~known).sum()) == terminal['n_pool'], 'Unresolved pool mismatch')
        require(int(known.sum()) == terminal['n_known'], 'Known marker pool mismatch')
        if terminal['dl_status'] == 'no_op_all_initially_known':
            require(terminal['n_pool'] == 0, 'Invalid all-known no-op')
        if terminal['dl_status'] == 'no_known_labels_archived_Undecided_terminal':
            require(terminal['n_known'] == 0, 'Invalid no-known-label endpoint')
        if terminal['dl_status'] == 'structural_insufficient_known_split':
            require(terminal['n_known'] == 1, 'Invalid insufficient-split endpoint')
        for stage, column, threshold in provider.STAGES:
            require(np.array_equal(calls.loc[known, 'initial'], calls.loc[known, column]),
                    'Known marker calls changed during refinement')
            scores = provider.score_predictions(calls[column], condition['library'], truth, scope)
            rows.append(dict(sample=task['sample'], patient=task['patient'], primary=task['primary'],
                budget=task['budget'], family=spec['family'], configuration=spec['configuration'],
                route=condition['route'], library=condition['library'], cutoff=str(condition['cutoff']),
                stage=stage, threshold=threshold, status='completed', terminal_valid=True,
                dl_status=terminal['dl_status'], training_executed=terminal['training_executed'],
                identical_result_reused=terminal['identical_result_reused'],
                n_cells=len(truth), n_training_classes=terminal['n_training_classes'],
                DL_features=terminal['DL_features'], terminal_manifest_sha256=terminal_hash,
                predictions_sha256=terminal['predictions_sha256'], truth_sha256=truth_hash,
                **scores))
        receipts[str(terminal_path / 'terminal_manifest.json')] = terminal_hash
        receipts[str(prediction_path)] = terminal['predictions_sha256']
    require(sha(spec_path) == spec_hash and sha(proof_path) == spec['parity_receipt']['sha256'],
            'Evaluation inputs changed while scoring')
    for path, digest in receipts.items():
        require(sha(path) == digest, f'Terminal artifacts changed while scoring: {path}')
    for relative, digest in lock['sources'].items():
        require(sha(provider.ROOT / relative) == digest, f'Evaluation source changed: {relative}')
    for path, digest in spec['source_hashes'].items():
        require(sha(path) == digest, f'Control source changed while scoring: {path}')
    require(sha(__file__) == evaluator_hash and sha(CORE) == CORE_SHA256 and
            sha(HERE / 'evaluate_core_lfine_sources.json') == LOCK_SHA256,
            'Evaluator changed while scoring')
    output.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(output / 'metrics.csv.gz', index=False)
    result = dict(status='completed', sample=task['sample'], budget=task['budget'],
        family=spec['family'], configuration=spec['configuration'],
        n_conditions=len(conditions), n_threshold_rows=len(rows), n_valid=len(rows),
        script_sha256=sha(__file__), specification_sha256=spec_hash,
        parity_receipt=spec['parity_receipt'], core_provider_sha256=CORE_SHA256,
        frozen_semantics=lock, truth_sha256=truth_hash, terminal_artifacts=receipts,
        outputs={'metrics.csv.gz': sha(output / 'metrics.csv.gz')},
        job=os.environ['SLURM_JOB_ID'], completed_at=datetime.now(timezone.utc).isoformat())
    (output / 'manifest.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    (output / 'COMPLETE').write_text(sha(output / 'manifest.json') + '\n')
    print(json.dumps({k: result[k] for k in ['status', 'family', 'n_conditions', 'n_threshold_rows']}))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spec', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    evaluate(args.spec, args.out)
