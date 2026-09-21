"""Fresh packaged-backend MLP controls on verified, immutable native-R inputs.

No truth or evaluation tables are read. Scientific refinement bodies are frozen
from the archived parameterized controls; only their common-module import moves.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import ast
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import traceback

HERE = Path(__file__).resolve().parent


def need(value, message):
    if not value:
        raise RuntimeError(message)


def load(path):
    return json.loads(Path(path).read_text())


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def write_new(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def check_receipt(directory, manifest, flag):
    need((directory / flag).read_text().strip() == sha(directory / manifest),
         f'Invalid receipt: {directory / flag}')


def check_sources(lock):
    need(lock['status'] == 'frozen_for_bounded_validation' and
         {'adapter', 'parameterized_refine', 'tasks', 'archived_refine', 'independent_verifier',
          'core_gate_validator', 'core_receipt_validator', 'extension_evaluator',
          'core_evaluator', 'core_evaluator_sources'}.issubset(lock['files']),
         'MLP adapter source lock is incomplete or not frozen')
    for name, record in lock['files'].items():
        need(sha(record['path']) == record['sha256'], f'Frozen adapter source changed: {name}')
    archived = Path(lock['files']['archived_refine']['path']).read_text()
    portable = (HERE / 'refine_controls_packaged.py').read_text()
    before = 'from ptc_common import BASE, RECOVERY, require_slurm, sha, utc, write_json, task_list, geometry_dir'
    after = 'from dgscrna.reference.backend.common import require_slurm, sha, utc, write_json'
    need(archived.count(before) == 1 and archived.replace(before, after) == portable,
         'Parameterized refinement differs beyond its common-module import')
    without_imports = lambda text: [ast.dump(n, include_attributes=False) for n in ast.parse(text).body
                                  if not isinstance(n, ast.ImportFrom)]
    need(without_imports(archived) == without_imports(portable), 'Scientific AST changed')


def verify_source(source, gate, gate_hash, task):
    campaign = Path(gate['output_root']).resolve()
    need(source.is_relative_to(campaign), 'Inputs must come from this new packaged campaign')
    check_receipt(source, 'run_manifest.json', 'COMPLETE')
    manifest = load(source / 'run_manifest.json')
    plan = load(source / 'run_config.json')
    proof = load(source / 'packaged_parity.json')
    need(manifest['status'] == 'completed' and manifest['terminal_condition_count'] == 192,
         'A full 192-condition packaged source is required')
    need(plan['config']['sample'] == task['sample'] and
         'hvg' + str(plan['config']['features']).removeprefix('hvg') == task['budget'],
         'Source sample or feature budget differs')
    need(plan['runtime'] == gate['runtime']['fingerprint'], 'Source runtime differs from the frozen wheel pilots')
    need(manifest['plan_sha256'] == sha(source / 'run_config.json'), 'Source configuration hash changed')
    for name, digest in plan['sources'].items():
        need(gate['installed_package']['files_sha256'].get('reference/' + name) == digest,
             'Source package differs from accepted wheel: ' + name)
    from dgscrna.reference.runner import _checked_receipt
    for stage in ['input', 'prepare', 'score', 'terminal']:
        need(_checked_receipt(source / 'checkpoints' / (stage + '.json'), source, manifest['plan_sha256']),
             'Source checkpoint invalid: ' + stage)
    need(proof['status'] == 'passed_exact' and proof['terminal_conditions'] == 192 and
         proof['run_manifest_sha256'] == sha(source / 'run_manifest.json') and
         proof['run_config_sha256'] == sha(source / 'run_config.json'), 'Source independent parity is invalid')
    if (source / 'PACKAGE_VERIFIED_COMPLETE').exists():
        check_receipt(source, 'package_acceptance.json', 'PACKAGE_VERIFIED_COMPLETE')
        accepted = load(source / 'package_acceptance.json')
        need(accepted['gate_sha256'] == gate_hash and accepted['parity_sha256'] == sha(source / 'packaged_parity.json'),
             'Core source belongs to a different release gate')
        manager = module('_frozen_core_receipts', HERE / 'manage_core.py')
        roster = load(gate['task_manifest']['path'])
        matching = next(t for t in roster if t['sample'] == task['sample'] and t['budget'] == task['budget'])
        need(source == campaign / 'core' / task['sample'] / task['budget'] and
             manager.acceptance(campaign, matching, gate_hash, gate['wheel']['sha256']) is not None,
             'Core source must include its full accepted Lfine evaluation')
    else:
        role = next((p for p in gate['pilots'] if p['role'] == 'full_roster'), None)
        need(role is not None and Path(role['path']).resolve().parent == source and
             role['sha256'] == sha(source / 'packaged_parity.json'), 'Unaccepted non-core source')
    need(manifest['reference_labels_used_for_fit'] is False, 'Source used truth for fitting')
    prep = Path(manifest['prepare']).resolve()
    need(prep == source / 'GBM' / task['sample'] / task['budget'], 'Unexpected source preparation path')
    snapshots = {str(source / n): sha(source / n) for n in
                 ['run_manifest.json', 'run_config.json', 'packaged_parity.json', 'COMPLETE']}
    pm = load(prep / 'prepare_manifest.json')
    binary = Path(pm['DL_binary']).resolve()
    need(binary.is_relative_to(source) and sha(binary) == pm['DL_binary_sha256'], 'Source DL matrix changed')
    for path in [binary, prep / 'cells.csv', prep / 'DL_features.txt', prep / 'prepare_manifest.json', prep / 'PREPARED']:
        snapshots[str(path)] = sha(path)
    for route in task['routes']:
        directory = prep / route
        check_receipt(directory, 'score_manifest.json', 'SCORE_COMPLETE')
        score = load(directory / 'score_manifest.json')
        need(Path(score['DL_binary']).resolve() == binary and score['DL_binary_sha256'] == sha(binary),
             'Scorer and preparation use different DL inputs')
        need(score['initial_sha256'] == sha(directory / 'initial_calls.csv.gz'), 'Initial labels changed')
        for name in ['score_manifest.json', 'SCORE_COMPLETE', 'cells.csv', 'initial_calls.csv.gz']:
            snapshots[str(directory / name)] = sha(directory / name)
    return prep, snapshots


def run(args):
    need(os.environ.get('SLURM_JOB_ID'), 'Scientific controls must run through SLURM')
    need(int(os.environ.get('SLURM_CPUS_PER_TASK', 0)) >= 4, 'Four CPUs are required for exact historical DL runtime')
    need(sys.flags.no_user_site and not sys.flags.optimize and not os.environ.get('PYTHONPATH'),
         'Use frozen installed Python -s, without PYTHONPATH or -O')
    os.chdir('/tmp')
    lock_path = HERE / 'mlp_controls_sources.json'
    lock = load(lock_path)
    check_sources(lock)
    worker = module('_frozen_core_gate', HERE / 'core_task.py')
    gate, _, package = worker.validate_gate(args.gate)
    tasks = load(HERE / 'mlp_controls_tasks.json')
    need(len(tasks) == 54 and 0 <= args.index < len(tasks), 'Invalid controls task roster/index')
    task = tasks[args.index]
    need(task['index'] == args.index and task['n_terminal_conditions'] == 8, 'Task identity mismatch')
    source = (args.source or (Path(gate['output_root']) / 'core' / task['sample'] / task['budget'])).resolve()
    prep, input_hashes = verify_source(source, gate, sha(args.gate), task)
    from dgscrna.reference.runner import doctor
    need(doctor(gate['runtime']['rscript'], gate['runtime'].get('reference_r_lib')) == gate['runtime']['fingerprint'],
         'Current runtime differs from accepted pilots')
    out = args.out.resolve()
    campaign = Path(gate['output_root']).resolve()
    need(out.is_relative_to(campaign / 'extensions') and not out.exists(),
         'Use a fresh output under campaign/extensions; existing artifacts must be preserved')
    os.environ['DGSCRNA_REFERENCE_OUT'] = str(out)
    os.environ['DGSCRNA_DL_CACHE_ROOT'] = str(out / 'DL_cache')
    os.environ['DGSCRNA_REQUIRE_SLURM'] = '1'
    out.mkdir(parents=True)
    started = dict(status='running', task=task, source=str(source), gate_sha256=sha(args.gate),
                   adapter_sources_sha256=sha(lock_path), package_root=str(package),
                   source_input_hashes=input_hashes, reference_labels_used_for_fit=False,
                   job=os.environ['SLURM_JOB_ID'], step=os.environ.get('SLURM_STEP_ID'),
                   started_at=datetime.now(timezone.utc).isoformat())
    write_new(out / 'run_config.json', started)
    try:
        from dgscrna.reference.backend import terminal
        terminal.threads()
        refine = module('_frozen_parameterized_refinement', HERE / 'refine_controls_packaged.py')
        refine.PARAMS.update(task['params'])
        terminal.refine = refine
        verifier = module('_independent_terminal_verifier', HERE / 'verify_packaged_parity.py')
        old = Path(gate['reference_root']) / 'GBM_DL_controls' / task['sample'] / task['budget'] / task['control']
        check_receipt(old, 'manifest.json', 'COMPLETE')
        om = load(old / 'manifest.json')
        need(om['params'] == task['params'] and om['n_conditions'] == 8 and
             om['refinement_source_sha256'] == lock['files']['archived_refine']['sha256'],
             'Archived control does not match this frozen scientific configuration')
        conditions, checks, caches, evaluation_conditions, pending_comparisons = [], [], set(), [], []
        for route in task['routes']:
            original, dest = prep / route, out / route
            dest.mkdir()
            for name in ['score_manifest.json', 'SCORE_COMPLETE', 'cells.csv', 'initial_calls.csv.gz']:
                shutil.copyfile(original / name, dest / name)
                need(sha(dest / name) == input_hashes[str(original / name)], 'Copied input differs')
            sm, em = load(dest / 'score_manifest.json'), load(old / route / 'score_manifest.json')
            actual_arms, expected_arms = verifier.arm_map(sm), verifier.arm_map(em)
            for library in task['libraries']:
                key = (library, task['cutoff'])
                aid, arm = actual_arms[key]
                eid, _ = expected_arms[key]
                terminal.finish_route(dest, only_arm=aid)
                td = dest / 'terminal' / aid
                tm = load(td / 'terminal_manifest.json')
                cache = Path(tm['cache_directory']).resolve()
                need(cache.is_relative_to(out / 'DL_cache'), 'Control borrowed an old/external DL cache')
                training = load(td / 'training_manifest.json')
                need(Path(training['provenance']['first_condition']).resolve().is_relative_to(out),
                     'Cache first condition belongs to an earlier run')
                caches.add(str(cache))
                record = dict(route=route, library=library, cutoff=task['cutoff'],
                              dl_status=tm['dl_status'], training_executed=tm['training_executed'])
                pending_comparisons.append((td, old / route / 'terminal' / eid, record, aid, eid))
                evaluation_conditions.append(dict(route=route, library=library, cutoff=task['cutoff'],
                                                  terminal_directory=str(td)))
        # Every new fit finishes before archived terminal arrays/models are opened.
        for actual, expected, record, aid, eid in pending_comparisons:
            conditions.append(verifier.compare_terminal(actual, expected, record, aid, eid, checks, fresh=False))
        need(len(conditions) == 8, 'Incomplete MLP controls roster')
        need(all(sha(path) == digest for path, digest in input_hashes.items()), 'Immutable source input changed')
        check_sources(lock)
        worker.validate_gate(args.gate)
        parity = dict(status='passed_exact', sample=task['sample'], budget=task['budget'],
                        control=task['control'], params=task['params'], terminal_conditions=8,
                        conditions=conditions, artifact_checks=checks, n_unique_fresh_caches=len(caches),
                        cache_paths_confined_to_unit=True, source_inputs_unchanged=True,
                        source_input_hashes=input_hashes, reference_labels_used_for_fit=False,
                        truth_opened_before_parity=False, L1_evaluation_opened=False, gate_sha256=sha(args.gate),
                        adapter_sources_sha256=sha(lock_path), config_sha256=sha(out / 'run_config.json'),
                        archived_control_manifest_sha256=sha(old / 'manifest.json'),
                        job=os.environ['SLURM_JOB_ID'], completed_at=datetime.now(timezone.utc).isoformat())
        write_new(out / 'parity.json', parity)
        spec = dict(sample=task['sample'], budget=task['budget'], family='MLP_only',
                    configuration=task['control'], expected_conditions=8, conditions=evaluation_conditions,
                    parity_receipt=dict(path=str(out / 'parity.json'), sha256=sha(out / 'parity.json')),
                    source_hashes={record['path']:record['sha256'] for record in lock['files'].values()})
        write_new(out / 'evaluation_spec.json', spec)
        command = [gate['runtime']['python'], '-s', str(HERE / 'evaluate_extension_lfine.py'),
                   '--spec', str(out / 'evaluation_spec.json'), '--out', str(out / 'evaluation')]
        with (out / 'evaluation.log').open('x') as log:
            result = subprocess.run(command, cwd='/tmp', stdout=log, stderr=subprocess.STDOUT)
        need(result.returncode == 0, 'Frozen Lfine evaluation failed; preserve output for review')
        check_receipt(out / 'evaluation', 'manifest.json', 'COMPLETE')
        evaluated = load(out / 'evaluation/manifest.json')
        need(evaluated['status'] == 'completed' and evaluated['n_conditions'] == 8 and
             evaluated['n_threshold_rows'] == evaluated['n_valid'] == 16 and
             evaluated['sample'] == task['sample'] and evaluated['budget'] == task['budget'] and
             evaluated['configuration'] == task['control'], 'Lfine result is incomplete or identifies another unit')
        need(evaluated['script_sha256'] == lock['files']['extension_evaluator']['sha256'] and
             evaluated['specification_sha256'] == sha(out / 'evaluation_spec.json') and
             evaluated['parity_receipt'] == spec['parity_receipt'] and
             evaluated['outputs']['metrics.csv.gz'] == sha(out / 'evaluation/metrics.csv.gz'),
             'Lfine result hashes differ from frozen adapter inputs')
        check_sources(lock)
        manifest = dict(parity, lfine_valid_threshold_rows=16, parity_sha256=sha(out / 'parity.json'),
                        lfine_manifest_sha256=sha(out / 'evaluation/manifest.json'),
                        lfine_metrics_sha256=sha(out / 'evaluation/metrics.csv.gz'))
        write_new(out / 'manifest.json', manifest)
        (out / 'COMPLETE').write_text(sha(out / 'manifest.json') + '\n')
        print(json.dumps(dict(status='passed_exact', output=str(out), terminal_conditions=8)), flush=True)
    except BaseException:
        write_new(out / 'FAILURE.json', dict(status='failed_preserved', traceback=traceback.format_exc(),
                                            job=os.environ['SLURM_JOB_ID']))
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gate', required=True, type=Path)
    parser.add_argument('--index', required=True, type=int)
    parser.add_argument('--source', type=Path, help='Accepted full packaged pilot; defaults to matching core unit')
    parser.add_argument('--out', required=True, type=Path)
    ARGS = parser.parse_args()
    ARGS.gate = ARGS.gate.resolve()
    run(ARGS)
