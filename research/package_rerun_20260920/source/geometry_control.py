"""Research adapter: native geometry with fresh, fixed-HVG2000 DL refinement.

Uses only newly accepted package core artifacts and the installed backend.
The old geometry terminal outputs are opened by the independent verifier only.
No old evaluation tables, labels, models or training caches enter fitting.
"""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import errno
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import uuid

import core_task


HERE = Path(__file__).resolve().parent
LIBRARIES = ('CM2_glioma_other', 'CM2_primary_all_context')
BUDGETS = ('hvg2000', 'hvg5000', 'all')
ROUTES = ('PCA30_SNN', 'PCA30_HDBSCAN_R', 'UMAP2_SNN', 'UMAP2_HDBSCAN_R')


def need(value, message):
    if not value:
        raise RuntimeError(message)


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def write_new(path, data):
    with Path(path).open('x') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')


def checked(path, manifest, flag):
    need((path / flag).read_text().strip() == sha(path / manifest), f'Invalid completion receipt: {path}')
    return load(path / manifest)


def file_inventory(root):
    """Detect core additions/deletions/rewrites without loading model objects."""
    return {str(p.relative_to(root)): (p.stat().st_size, p.stat().st_mtime_ns)
            for p in sorted(root.rglob('*')) if p.is_file()}


def source_run(path, gate, sample, budget, *, pilot=False):
    from dgscrna.reference.runner import _checked_receipt
    path = path.resolve()
    run = checked(path, 'run_manifest.json', 'COMPLETE')
    plan = load(path / 'run_config.json')
    need(run['status'] == 'completed' and run['terminal_condition_count'] == 192,
         'Geometry controls require the complete fresh 192-arm core source')
    need(run['plan_sha256'] == sha(path / 'run_config.json'), 'Source run plan hash changed')
    config = plan['config']
    need(config['sample'] == sample and ('all' if config['features'] == 'all' else 'hvg' + config['features']) == budget,
         'Source sample/budget differs from selected control')
    need(config['route'] == config['library'] == config['cutoff'] == 'all', 'Source has a reduced condition roster')
    need(plan['runtime'] == gate['runtime']['fingerprint'], 'Source runtime differs from release pilots')
    for key, digest in plan['sources'].items():
        need(gate['installed_package']['files_sha256'].get('reference/' + key) == digest,
             f'Source package differs: {key}')
    for stage in ('input', 'prepare', 'score', 'terminal'):
        need(_checked_receipt(path / 'checkpoints' / (stage + '.json'), path, run['plan_sha256']),
             f'Missing source checkpoint: {stage}')
    parity = load(path / 'packaged_parity.json')
    need(parity['status'] == 'passed_exact' and parity['terminal_conditions'] == 192
         and parity['sample'] == sample and parity['budget'] == budget,
         'Source does not have full independent exact parity')
    need(parity['run_manifest_sha256'] == sha(path / 'run_manifest.json')
         and parity['run_config_sha256'] == sha(path / 'run_config.json'), 'Source parity metadata changed')
    if pilot:
        record = next(row for row in gate['pilots'] if row['role'] == 'full_roster')
        need(Path(record['path']).resolve() == path / 'packaged_parity.json'
             and record['sha256'] == sha(path / 'packaged_parity.json'), 'Pilot source is not the frozen full-roster pilot')
        authorization = dict(kind='frozen_full_roster_pilot', proof_sha256=record['sha256'])
    else:
        acceptance = checked(path, 'package_acceptance.json', 'PACKAGE_VERIFIED_COMPLETE')
        need(acceptance['status'] == 'package_run_independently_verified'
             and acceptance['gate_sha256'] == gate['_gate_sha256']
             and acceptance['run_manifest_sha256'] == sha(path / 'run_manifest.json')
             and acceptance['parity_sha256'] == sha(path / 'packaged_parity.json'),
             'Source package acceptance differs from the published gate')
        from manage_core import acceptance as verify_core_acceptance
        tasks = load(core_task.check_file(gate['task_manifest']))
        task = next(t for t in tasks if t['sample'] == sample and t['budget'] == budget)
        need(verify_core_acceptance(Path(gate['output_root']).resolve(), task,
                                    gate['_gate_sha256'], gate['wheel']['sha256']) is not None,
             'Source core lacks its verified final Lfine evaluation')
        authorization = dict(kind='accepted_new_package_core', acceptance_sha256=sha(path / 'package_acceptance.json'))
    prep = Path(run['prepare']).resolve()
    need(prep == path / 'GBM' / sample / budget, 'Unexpected source preparation path')
    checked(prep, 'prepare_manifest.json', 'PREPARED')
    return prep, dict(root=str(path), sample=sample, budget=budget,
                      run_manifest_sha256=sha(path / 'run_manifest.json'),
                      run_config_sha256=sha(path / 'run_config.json'),
                      parity_sha256=sha(path / 'packaged_parity.json'), **authorization)


def arm_map(manifest):
    mapping = {}
    for aid, arm in manifest['arms'].items():
        key = (arm['library'], str(arm['cutoff']))
        need(key not in mapping, f'Duplicate marker arm: {key}')
        mapping[key] = aid
    return mapping


def copy_input(source, destination, recorded, *, binary_link=False):
    need(source.is_file(), f'Missing immutable input: {source}')
    digest = sha(source)
    mode = 'copy'
    if binary_link:
        try:
            os.link(source, destination)
            mode = 'read_only_hardlink'
        except OSError as error:
            if error.errno != errno.EXDEV:
                raise
            shutil.copyfile(source, destination)
    else:
        shutil.copyfile(source, destination)
    need(sha(destination) == digest, f'Input copy differs: {source}')
    recorded.append(dict(source=str(source), destination=str(destination), sha256=digest, mode=mode))


def load_verifier(gate):
    path = core_task.check_file(gate['launcher_files']['verify_packaged_parity.py'])
    spec = importlib.util.spec_from_file_location('geometry_independent_terminal_verifier', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gate', required=True, type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--index', type=int, help='0..362 in the frozen core roster filtered to HVG2000/HVG5000/all')
    mode.add_argument('--pilot-full-run', type=Path, help='Explicit frozen TKU4163 full-roster HVG2000 pilot source')
    args = parser.parse_args()
    need(os.environ.get('SLURM_JOB_ID'), 'Scientific work requires SLURM')
    need(sys.flags.no_user_site and not sys.flags.optimize and not os.environ.get('PYTHONPATH'),
         'Use the frozen installed Python with -s, without -O or PYTHONPATH')
    need(int(os.environ.get('SLURM_CPUS_PER_TASK', 0)) >= 4, 'Geometry parity requires the original four Torch threads')
    os.chdir('/tmp')
    gate_path = args.gate.resolve()
    gate, tasks, package = core_task.validate_gate(gate_path)
    from dgscrna.reference.runner import doctor
    need(doctor(gate['runtime']['rscript'], gate['runtime'].get('reference_r_lib')) == gate['runtime']['fingerprint'],
         'Current installed runtime differs from the published release pilots')
    gate_hash = sha(gate_path)
    gate['_gate_sha256'] = gate_hash
    root = Path(gate['output_root']).resolve()
    geometry_tasks = [t for t in tasks if t['budget'] in BUDGETS]
    need(len(geometry_tasks) == 363, 'Geometry control roster must contain 363 sample/budget units')
    pilot = args.pilot_full_run is not None
    if pilot:
        sample, budget = 'TKU4163', 'hvg2000'
        source_root = fixed_root = args.pilot_full_run.resolve()
        dest = root / 'pilots/geometry_fixed_DL2000_TKU4163_hvg2000'
        task_index = None
    else:
        need(0 <= args.index < len(geometry_tasks), 'Invalid geometry task index')
        task = geometry_tasks[args.index]
        sample, budget, task_index = task['sample'], task['budget'], task['index']
        source_root = root / 'core' / sample / budget
        fixed_root = root / 'core' / sample / 'hvg2000'
        dest = root / 'extensions/geometry_fixed_DL2000' / sample / budget
    dest = dest.resolve()
    need(not dest.exists(), f'Preserve the existing extension attempt and inspect it before retry: {dest}')
    need(dest.is_relative_to(root) and not dest.is_relative_to((root / 'core').resolve())
         and not dest.is_relative_to(source_root.resolve()) and not dest.is_relative_to(fixed_root.resolve()),
         'Extension must stay in its campaign and must not modify a core/source directory')
    geometry, geometry_record = source_run(source_root, gate, sample, budget, pilot=pilot)
    if fixed_root.resolve() == source_root.resolve():
        fixed, fixed_record = geometry, dict(geometry_record)
    else:
        fixed, fixed_record = source_run(fixed_root, gate, sample, 'hvg2000', pilot=pilot)
    fixed_manifest = checked(fixed, 'prepare_manifest.json', 'PREPARED')
    need(fixed_manifest['features']['DL'] == 2000 and fixed_manifest['assay'] == 'RNA',
         'Fixed DL input is not the native normalized RNA HVG2000 matrix')
    original = Path(gate['reference_root']).resolve() / 'GBM' / sample / budget
    old_fixed = Path(gate['reference_root']).resolve() / 'GBM' / sample / 'hvg2000'
    need(sha(fixed / 'DL.float32.bin') == fixed_manifest['DL_binary_sha256']
         == sha(old_fixed / 'DL.float32.bin'), 'Fixed new/reference DL matrices differ')
    need(sha(fixed / 'cells.csv') == sha(geometry / 'cells.csv') == sha(old_fixed / 'cells.csv'),
         'Geometry and fixed-DL cell order differ')
    source_inventory = {str(p): file_inventory(p) for p in {source_root, fixed_root}}
    selected = []
    for route in ROUTES:
        source_manifest = checked(geometry / route, 'score_manifest.json', 'SCORE_COMPLETE')
        old_manifest = checked(original / route, 'score_manifest.json', 'SCORE_COMPLETE')
        actual_arms, expected_arms = arm_map(source_manifest), arm_map(old_manifest)
        for library in LIBRARIES:
            key = (library, 'mean')
            need(key in actual_arms and key in expected_arms, f'Missing prespecified geometry arm: {route}/{library}')
            expected = original / route / 'terminal_geometry_only_DL2000' / expected_arms[key]
            checked(expected, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
            selected.append(dict(route=route, library=library, cutoff='mean', arm=actual_arms[key],
                                 reference_arm=expected_arms[key], reference_terminal=str(expected)))
    need(len(selected) == 8, 'Geometry comparison requires exactly eight terminal conditions')
    dest.mkdir(parents=True, exist_ok=False)
    source_hash = sha(__file__)
    evaluator = HERE / 'evaluate_extension_lfine.py'
    evaluator_hash = sha(evaluator)
    protocol = dict(status='fresh_geometry_control_started', sample=sample, budget=budget,
                    geometry_index=args.index, core_task_index=task_index, pilot_only=pilot,
                    expected_full_units=363, gate_sha256=gate_hash, adapter_sha256=source_hash,
                    extension_evaluator_sha256=evaluator_hash,
                    installed_package_root=str(package), wheel_sha256=gate['wheel']['sha256'],
                    source_geometry=geometry_record, source_fixed_DL=fixed_record,
                    fixed_DL_features=2000, fixed_DL_sha256=fixed_manifest['DL_binary_sha256'],
                    libraries=list(LIBRARIES), cutoff='mean', routes=list(ROUTES),
                    scientific_contract='Use each geometry-specific native RNA scoring result; replace only DL expression input with the same sample normalized RNA HVG2000 matrix',
                    no_new_DEG_calculation=True, no_reference_truth_for_fit=True,
                    old_L1_evaluation_used=False, old_training_cache_used=False,
                    started_at=utc(), job=os.environ['SLURM_JOB_ID'])
    write_new(dest / 'protocol.json', protocol)
    try:
        copied = []
        local_fixed = dest / 'inputs/fixed_hvg2000'
        local_fixed.mkdir(parents=True)
        for filename in ['cells.csv', 'DL_features.txt']:
            copy_input(fixed / filename, local_fixed / filename, copied)
        copy_input(fixed / 'DL.float32.bin', local_fixed / 'DL.float32.bin', copied, binary_link=True)
        local_manifest = dict(fixed_manifest, DL_binary=str(local_fixed / 'DL.float32.bin'))
        write_new(local_fixed / 'prepare_manifest.json', local_manifest)
        (local_fixed / 'PREPARED').write_text(sha(local_fixed / 'prepare_manifest.json') + '\n')
        for route in ROUTES:
            local_source = dest / 'scores' / route
            local_source.mkdir(parents=True)
            for filename in ['score_manifest.json', 'SCORE_COMPLETE', 'initial_calls.csv.gz', 'cells.csv', 'clusters.csv']:
                copy_input(geometry / route / filename, local_source / filename, copied)
        write_new(dest / 'input_provenance.json', dict(files=copied, local_fixed_prepare_sha256=sha(local_fixed / 'prepare_manifest.json'),
                                                    prepared_manifest_change='Only DL_binary path redirected to the exact copied/linked input'))
        environment = os.environ.copy()
        for name in ['PYTHONPATH', 'DGSCRNA_DL_CACHE_ROOT', 'DGSCRNA_SITE_ROOT', 'DGSCRNA_REFERENCE_OUT', 'DGSCRNA_EXAMPLE_OUT']:
            environment.pop(name, None)
        environment.update(DGSCRNA_REFERENCE_OUT=str(dest), DGSCRNA_REQUIRE_SLURM='1',
                           DGSCRNA_DL_CACHE_ROOT=str(dest / 'DL_cache' / uuid.uuid4().hex),
                           PYTHONNOUSERSITE='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
        logs = dest / 'logs'
        logs.mkdir()
        for record in selected:
            need(sha(gate_path) == gate_hash and sha(__file__) == source_hash, 'Frozen gate/adapter changed during execution')
            command = [gate['runtime']['python'], '-s', '-m', 'dgscrna.reference.backend.terminal',
                       str(dest / 'scores' / record['route']), record['arm'], str(local_fixed)]
            with (logs / (record['route'] + '_' + record['arm'] + '.log')).open('x') as stream:
                result = subprocess.run(command, env=environment, cwd='/tmp', stdout=stream, stderr=subprocess.STDOUT)
            need(result.returncode == 0, f'Fresh geometry refinement failed: {record}')
        # Verification opens archived terminal arrays/weights only after every
        # fresh fit has completed. These reference artifacts cannot enter fit.
        verifier = load_verifier(gate)
        conditions, checks = [], []
        for record in selected:
            actual = dest / 'scores' / record['route'] / 'terminal_geometry_only_DL2000' / record['arm']
            terminal = checked(actual, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
            need(terminal['family'] == 'geometry_only_fixed_DL2000', 'Wrong terminal family')
            cache = Path(terminal['cache_directory']).resolve()
            need(cache.is_relative_to(dest / 'DL_cache'), 'A terminal cache escaped this fresh extension unit')
            training = load(actual / 'training_manifest.json')
            need(Path(training['provenance']['first_condition']).resolve().is_relative_to(dest),
                 'Training reused a cache originating outside this fresh unit')
            condition = dict(record, dl_status=terminal['dl_status'], training_executed=terminal['training_executed'])
            conditions.append(verifier.compare_terminal(actual, Path(record['reference_terminal']), condition,
                                                        record['arm'], record['reference_arm'], checks, fresh=False))
        for source, snapshot in source_inventory.items():
            need(file_inventory(Path(source)) == snapshot, f'Accepted core files changed: {source}')
        # Recheck actual immutable inputs after training, including the shared binary.
        for item in copied:
            need(sha(item['source']) == item['sha256'] == sha(item['destination']), 'A fitting input changed during execution')
        accepted = dict(status='passed_exact', family='geometry_only_fixed_DL2000', sample=sample, budget=budget,
                        terminal_conditions=8, trained_conditions=sum(c['training_executed'] for c in conditions),
                        valid_no_training_conditions=sum(not c['training_executed'] for c in conditions),
                        gate_sha256=gate_hash, adapter_sha256=source_hash, protocol_sha256=sha(dest / 'protocol.json'),
                        input_provenance_sha256=sha(dest / 'input_provenance.json'),
                        conditions=conditions, checks=checks, numerical_tolerance=dict(rtol=0, atol=0),
                        core_directories_unchanged=True, caches_confined_to_fresh_extension_unit=True,
                        no_reference_truth_read=True, no_old_L1_evaluation=True,
                        pilot_only=pilot, job=os.environ['SLURM_JOB_ID'], completed_at=utc())
        write_new(dest / 'geometry_parity.json', accepted)
        need(sha(evaluator) == evaluator_hash, 'Extension evaluator changed during execution')
        evaluation_conditions = [dict(route=r['route'], library=r['library'], cutoff=r['cutoff'],
                                      terminal_directory=str(dest / 'scores' / r['route'] /
                                                             'terminal_geometry_only_DL2000' / r['arm']))
                                 for r in selected]
        evaluation_spec = dict(sample=sample, budget=budget, family='geometry_only_fixed_DL2000',
                               configuration='fixed_normalized_RNA_HVG2000_DL', expected_conditions=8,
                               conditions=evaluation_conditions,
                               parity_receipt=dict(path=str(dest / 'geometry_parity.json'),
                                                   sha256=sha(dest / 'geometry_parity.json')),
                               source_hashes={str(Path(__file__).resolve()): source_hash,
                                              str(evaluator): evaluator_hash,
                                              gate['launcher_files']['verify_packaged_parity.py']['path']:
                                                  gate['launcher_files']['verify_packaged_parity.py']['sha256'],
                                              gate['launcher_files']['evaluate_core_lfine.py']['path']:
                                                  gate['launcher_files']['evaluate_core_lfine.py']['sha256'],
                                              gate['launcher_files']['evaluate_core_lfine_sources.json']['path']:
                                                  gate['launcher_files']['evaluate_core_lfine_sources.json']['sha256']})
        write_new(dest / 'evaluation_spec.json', evaluation_spec)
        with (logs / 'Lfine_evaluation.log').open('x') as stream:
            result = subprocess.run([gate['runtime']['python'], '-s', str(evaluator),
                                     '--spec', str(dest / 'evaluation_spec.json'), '--out', str(dest / 'evaluation')],
                                    env=environment, cwd='/tmp', stdout=stream, stderr=subprocess.STDOUT)
        need(result.returncode == 0, 'Lfine control evaluation failed; preserve the parity proof and inspect its log')
        evaluation = checked(dest / 'evaluation', 'manifest.json', 'COMPLETE')
        need(evaluation['status'] == 'completed' and evaluation['n_conditions'] == 8
             and evaluation['n_threshold_rows'] == evaluation['n_valid'] == 16,
             'Geometry Lfine evaluation is incomplete')
        need(sha(dest / 'evaluation/metrics.csv.gz') == evaluation['outputs']['metrics.csv.gz'],
             'Lfine metric payload changed')
        for source, snapshot in source_inventory.items():
            need(file_inventory(Path(source)) == snapshot, f'Core changed during Lfine evaluation: {source}')
        accepted.update(parity_receipt=evaluation_spec['parity_receipt'],
                        lfine_manifest=str(dest / 'evaluation/manifest.json'),
                        lfine_manifest_sha256=sha(dest / 'evaluation/manifest.json'),
                        lfine_metrics_sha256=evaluation['outputs']['metrics.csv.gz'], lfine_valid_threshold_rows=16)
        write_new(dest / 'geometry_acceptance.json', accepted)
        (dest / 'GEOMETRY_VERIFIED_COMPLETE').write_text(sha(dest / 'geometry_acceptance.json') + '\n')
        print(json.dumps({k: accepted[k] for k in ['status', 'sample', 'budget', 'terminal_conditions', 'trained_conditions', 'pilot_only']}), flush=True)
    except BaseException as error:
        write_new(dest / 'failure_preserved.json', dict(status='failed_preserved', error=repr(error),
                                                       gate_sha256=gate_hash, adapter_sha256=source_hash,
                                                       job=os.environ['SLURM_JOB_ID'], at=utc()))
        raise


if __name__ == '__main__':
    main()
