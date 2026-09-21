"""Fresh packaged A2 seed-mechanism controls, followed by exact parity and Lfine."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
sys.path.insert(0, str(PARENT))
import core_task
import geometry_control as shared


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def call(command, environment, logfile):
    with logfile.open('x') as stream:
        result = subprocess.run(command, cwd='/tmp', env=environment,
                                stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'Stage failed ({result.returncode}); preserved log: {logfile}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gate', required=True, type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--index', type=int)
    mode.add_argument('--pilot-full-run', type=Path)
    args = parser.parse_args()
    assert os.environ.get('SLURM_JOB_ID') and int(os.environ.get('SLURM_CPUS_PER_TASK', 0)) >= 4
    assert sys.flags.no_user_site and not sys.flags.optimize and not os.environ.get('PYTHONPATH')
    os.chdir('/tmp')
    lock = load(HERE / 'SOURCE_LOCK.json')
    assert lock['status'] == 'reviewed_ready_for_pilot'
    for path, digest in lock['files'].items():
        assert sha(path) == digest, path
    gate, tasks, package = core_task.validate_gate(args.gate.resolve())
    gate_hash = sha(args.gate)
    gate['_gate_sha256'] = gate_hash
    from dgscrna.reference.runner import doctor
    assert doctor(gate['runtime']['rscript'], gate['runtime'].get('reference_r_lib')) == gate['runtime']['fingerprint']
    root = Path(gate['output_root'])
    tasks = [t for t in tasks if t['budget'] in {'hvg2000', 'hvg5000'}]
    assert len(tasks) == 242
    pilot = args.pilot_full_run is not None
    if pilot:
        sample, budget = 'TKU4163', 'hvg2000'
        source = args.pilot_full_run.resolve()
        output = root / 'pilots/no_cluster_TKU4163_hvg2000'
    else:
        assert 0 <= args.index < len(tasks)
        sample, budget = tasks[args.index]['sample'], tasks[args.index]['budget']
        source = root / 'core' / sample / budget
        output = root / 'extensions/no_cluster' / sample / budget
    output = output.resolve()
    assert output.is_relative_to(root.resolve()) and not output.is_relative_to((root / 'core').resolve())
    assert not output.exists(), f'Preserve existing attempt: {output}'
    prep, authorization = shared.source_run(source, gate, sample, budget, pilot=pilot)
    before = shared.file_inventory(source)
    pm = shared.checked(prep, 'prepare_manifest.json', 'PREPARED')
    assert pm['assay'] == 'RNA'
    original_protocol = load(HERE / 'protocol_original.json')
    derivation = load(HERE / 'SOURCE_DERIVATION.json')
    assert sha(Path(derivation['original_protocol'])) == derivation['original_protocol_sha256']
    assert original_protocol == load(Path(derivation['original_protocol']))
    assert sample in original_protocol['samples'] and budget in original_protocol['budgets']
    arms = original_protocol['arms']
    assert [a['lambda'] for a in arms] == [0, .5, 1, 1.5, 2]
    assert sha(Path(gate['markers']['path'])) == original_protocol['marker_sha256']
    output.mkdir(parents=True)
    route = output / 'cellwise_seed'
    configuration = dict(sample=sample, budget=budget, prep=str(prep), dest=str(output), route_dir=str(route),
        marker_file=gate['markers']['path'], marker_sha256=gate['markers']['sha256'],
        protocol_sha256=derivation['original_protocol_sha256'], source_bundle_sha256=sha(HERE / 'SOURCE_LOCK.json'),
        input=dict(prepare_manifest_sha256=sha(prep / 'prepare_manifest.json'),
                   expression_sha256=sha(prep / 'expression_PCA30.rds'), cells_sha256=sha(prep / 'cells.csv')),
        arms=arms)
    configuration['input_signature'] = hashlib.sha256(json.dumps(configuration, sort_keys=True).encode()).hexdigest()
    write(output / 'config.json', configuration)
    write(output / 'run_protocol.json', dict(status='started', sample=sample, budget=budget,
        family='no_cluster_seed_control', gate_sha256=gate_hash, wheel_sha256=gate['wheel']['sha256'],
        source=authorization, source_lock_sha256=sha(HERE / 'SOURCE_LOCK.json'),
        seed_mechanism_replaced=True, old_training_cache_used=False, reference_labels_used_for_fit=False,
        pilot=pilot, expected_units=242, terminal_conditions=5, job=os.environ['SLURM_JOB_ID']))
    environment = os.environ.copy()
    for name in ['PYTHONPATH', 'PYTHONHOME', 'R_HOME', 'R_LIBS', 'R_LIBS_USER', 'R_LIBS_SITE',
                 'DGSCRNA_REFERENCE_R_LIB', 'DGSCRNA_DL_CACHE_ROOT', 'DGSCRNA_EXAMPLE_OUT']:
        environment.pop(name, None)
    environment.update(DGSCRNA_REFERENCE_OUT=str(output), DGSCRNA_DL_CACHE_ROOT=str(output / 'DL_cache'),
        DGSCRNA_REQUIRE_SLURM='1', R_ENVIRON_USER='/dev/null', R_PROFILE_USER='/dev/null',
        R_LIBS_USER='', R_LIBS_SITE='', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    try:
        call([gate['runtime']['rscript'], '--vanilla', str(HERE / 'score_cellwise.R'), str(output / 'config.json')],
             environment, output / 'score.log')
        score = shared.checked(route, 'score_manifest.json', 'SCORE_COMPLETE')
        assert set(score['arms']) == {a['id'] for a in arms}
        assert score['no_cluster_or_DEG_used'] and score['seed_mechanism_replaced']
        call([gate['runtime']['python'], '-s', '-m', 'dgscrna.reference.backend.terminal', str(route), 'all'],
             environment, output / 'terminal.log')
        # Open prior fitted outputs only after all five fresh endpoints exist.
        original = Path(derivation['original_protocol']).parent / 'GBM' / sample / budget / 'cellwise_seed'
        verifier = shared.load_verifier(gate)
        checks, comparisons = [], []
        original_score = shared.checked(original, 'score_manifest.json', 'SCORE_COMPLETE')
        assert score['arms'] == original_score['arms']
        for name in ['cells.csv', 'initial_calls.csv.gz', 'cellwise_diagnostics.csv',
                     'seed_arm_counts.csv', 'marker_retention.csv']:
            verifier.frame_equal(verifier.frame(route / name), verifier.frame(original / name), name, checks)
        call([gate['runtime']['rscript'], '--vanilla', str(HERE / 'verify_scores.R'),
              str(route / 'cellwise_scores.rds'), str(original / 'cellwise_scores.rds')],
             environment, output / 'R_parity.log')
        for aid, arm in score['arms'].items():
            terminal = route / 'terminal' / aid
            tm = shared.checked(terminal, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
            cache = Path(tm['cache_directory']).resolve()
            assert cache.is_relative_to(output / 'DL_cache')
            train = load(terminal / 'training_manifest.json')
            assert Path(train['provenance']['first_condition']).resolve().is_relative_to(output)
            record = dict(route='cellwise_seed', library=arm['library'], cutoff=str(arm['cutoff']),
                          dl_status=tm['dl_status'], training_executed=tm['training_executed'])
            comparisons.append(verifier.compare_terminal(terminal, original / 'terminal' / aid,
                record, aid, aid, checks, fresh=False))
        parity = dict(status='passed_exact', sample=sample, budget=budget, conditions=comparisons,
            checks=checks, cellwise_scores_exact=True, original_lambda_grid_preserved=True,
            terminal_conditions=5, original_protocol_sha256=derivation['original_protocol_sha256'],
            source_lock_sha256=sha(HERE / 'SOURCE_LOCK.json'), job=os.environ['SLURM_JOB_ID'])
        write(output / 'parity.json', parity)
        conditions = [dict(route='cellwise_seed', library=arm['library'], cutoff=str(arm['cutoff']),
                           terminal_directory=str(route / 'terminal' / aid)) for aid, arm in score['arms'].items()]
        spec = dict(sample=sample, budget=budget, family='no_cluster_seed_control', configuration='lambda_grid',
            expected_conditions=5, conditions=conditions,
            parity_receipt=dict(path=str(output / 'parity.json'), sha256=sha(output / 'parity.json')),
            source_hashes=lock['files'])
        write(output / 'evaluation_spec.json', spec)
        call([gate['runtime']['python'], '-s', str(PARENT / 'evaluate_extension_lfine.py'),
              '--spec', str(output / 'evaluation_spec.json'), '--out', str(output / 'evaluation')],
             environment, output / 'evaluation.log')
        evaluation = shared.checked(output / 'evaluation', 'manifest.json', 'COMPLETE')
        assert evaluation['n_conditions'] == 5 and evaluation['n_valid'] == 10
        assert shared.file_inventory(source) == before, 'Source core changed'
        for path, digest in lock['files'].items():
            assert sha(path) == digest, path
        accepted = dict(status='passed_exact', sample=sample, budget=budget, gate_sha256=gate_hash,
            terminal_conditions=5, lfine_rows=10, source_core_unchanged=True,
            source_lock_sha256=sha(HERE / 'SOURCE_LOCK.json'), parity_sha256=sha(output / 'parity.json'),
            evaluation_manifest_sha256=sha(output / 'evaluation/manifest.json'),
            job=os.environ['SLURM_JOB_ID'], completed_at=datetime.now(timezone.utc).isoformat())
        write(output / 'acceptance.json', accepted)
        (output / 'COMPLETE').write_text(sha(output / 'acceptance.json') + '\n')
        print(json.dumps(accepted))
    except BaseException as error:
        write(output / 'failure.json', dict(error=repr(error), artifacts_preserved=True,
              job=os.environ['SLURM_JOB_ID'], time=datetime.now(timezone.utc).isoformat()))
        raise


if __name__ == '__main__':
    main()
