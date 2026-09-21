"""SLURM worker: frozen dependencies -> fresh adapter -> hash-bound acceptance."""
from pathlib import Path
import argparse
import os
import socket
import subprocess
import sys
import traceback
import common as c
import receipts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', required=True, type=Path)
    parser.add_argument('--index', type=int)
    args = parser.parse_args()
    c.safe_identity()
    c.need(sys.flags.no_user_site and not os.environ.get('PYTHONPATH')
           and int(os.environ.get('SLURM_CPUS_PER_TASK', 0)) >= 4,
           'Use installed Python -s with four SLURM CPUs and no PYTHONPATH')
    gate_path = args.campaign.resolve()
    gate_hash = c.sha(gate_path)
    gate, tasks, release, root = c.gate_metadata(gate_path)
    index = args.index if args.index is not None else int(os.environ['SLURM_ARRAY_TASK_ID'])
    c.need(0 <= index < len(tasks), 'Invalid extension task index')
    task = tasks[index]
    unit = c.worker_receipt_root(root, index)
    # A directory is a durable no-retry claim. A failed or interrupted worker
    # requires explicit review; neither a scheduler requeue nor manual duplicate
    # may launch another fit into its outputs.
    unit.parent.mkdir(parents=True, exist_ok=True)
    unit.mkdir()
    attempt = dict(status='preflight', task=task, campaign_gate_sha256=gate_hash,
                   release_gate_sha256=gate['release_gate']['sha256'],
                   worker_sha256=c.sha(__file__), started_at=c.utc(), pid=os.getpid(), host=socket.gethostname(),
                   job=os.environ['SLURM_JOB_ID'], array_job=os.environ.get('SLURM_ARRAY_JOB_ID'),
                   array_task=os.environ.get('SLURM_ARRAY_TASK_ID'), scientific_stage_entered=False)
    c.write(unit / 'attempt.json', attempt)
    try:
        c.need(str(attempt['array_task']) == str(index) and attempt['array_job'], 'Registered array worker identity required')
        state = c.load(root / 'control/extension_manager/state.json')
        c.need(state['campaign_gate_sha256'] == gate_hash, 'Controller state belongs to another campaign')
        registered = any(s['status'] == 'submitted' and s['job_id'] == attempt['array_job'] and index in s['indices']
                         for s in state['submissions'])
        intent = state.get('submission_intent', {})
        c.need(registered or (intent.get('name') == os.environ.get('SLURM_JOB_NAME') and index in intent.get('indices', [])),
               'Worker is not bound to a durable registered submission/intent')
        dependencies = c.dependency_proofs(task, release, root)
        c.need(dependencies is not None, 'Core dependencies are not PACKAGE_VERIFIED_COMPLETE')
        output = Path(task['output']).resolve()
        c.need(output.is_relative_to(root / 'extensions') and not output.exists(),
               f'Preserve existing or unsafe scientific output: {output}')
        for dependency in task['extension_dependencies']:
            unit_dependency = c.worker_receipt_root(root, dependency)
            accepted_dependency = c.receipt(unit_dependency, 'acceptance.json', 'EXTENSION_VERIFIED_COMPLETE')
            c.need(accepted_dependency['status'] == 'extension_independently_verified'
                   and accepted_dependency['campaign_gate_sha256'] == gate_hash
                   and accepted_dependency['task_index'] == dependency
                   and accepted_dependency['sample'] == task['sample'] and accepted_dependency['budget'] == task['budget']
                   and tasks[dependency]['configuration']['space'] == 'noDR', 'A1 noDR dependency differs')
            c.need(accepted_dependency['artifact_stats'] == c.artifact_stats(accepted_dependency['result']['artifacts']),
                   'Accepted A1 noDR artifacts changed before sibling geometry use')
        old_reference = c.historical_A1_ready(task, deep=True)
        c.need(old_reference is not None, 'Original A1 per-space reference is incomplete')
        environment = os.environ.copy()
        for name in ['PYTHONPATH', 'PYTHONHOME', 'R_HOME', 'R_LIBS', 'R_LIBS_USER', 'R_LIBS_SITE',
                     'DGSCRNA_REFERENCE_OUT', 'DGSCRNA_DL_CACHE_ROOT', 'DGSCRNA_ONLY_ROUTE']:
            environment.pop(name, None)
        environment.update(PYTHONNOUSERSITE='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                           MKL_NUM_THREADS='1', NUMBA_NUM_THREADS='1', R_ENVIRON_USER='/dev/null', R_PROFILE_USER='/dev/null')
        if task['family'] == 'A1':
            runtime = gate['A1_runtime']['python']
            command = [runtime, '-s', task['script'], '--package-run', task['configuration']['package_run'],
                       '--out', str(root / 'extensions/A1'), '--rscript', release['runtime']['rscript'],
                       '--spaces', task['configuration']['space']]
            c.need(not Path(task['verification_output']).exists(), 'Preserve existing A1 independent verification')
        else:
            command = [release['runtime']['python'], '-s', task['script'], '--gate', gate['release_gate']['path'], *task['arguments']]
        attempt.update(status='scientific_started', scientific_stage_entered=True, dependencies=dependencies,
                       command=command, original_reference=old_reference, scientific_started_at=c.utc())
        c.write(unit / 'attempt.json', attempt)
        with (unit / 'adapter.log').open('x') as stream:
            result = subprocess.run(command, env=environment, cwd='/tmp', stdout=stream, stderr=subprocess.STDOUT)
        c.need(result.returncode == 0, f'Adapter failed ({result.returncode}); inspect {unit / "adapter.log"}')
        if task['family'] == 'A1':
            verify = [gate['A1_runtime']['python'], '-s', str(c.PARENT / 'embedding_adapter/verify.py'),
                      '--actual', str(output.parent), '--expected', str(Path(task['configuration']['expected_reference']).parent),
                      '--rscript', release['runtime']['rscript'], '--evaluation-python', release['runtime']['python'],
                      '--out', task['verification_output'], '--spaces', task['configuration']['space']]
            attempt.update(status='independent_verification_started', verification_command=verify)
            c.write(unit / 'attempt.json', attempt)
            with (unit / 'independent_verification.log').open('x') as stream:
                result = subprocess.run(verify, env=environment, cwd='/tmp', stdout=stream, stderr=subprocess.STDOUT)
            c.need(result.returncode == 0, 'A1 independent verification failed; preserve fresh fit and logs')
        proof = receipts.validate(task, gate['release_gate']['sha256'], deep=True)
        c.need(c.sha(gate_path) == gate_hash, 'Campaign gate changed in flight')
        c.gate_metadata(gate_path)
        c.need(c.dependency_proofs(task, release, root) == dependencies, 'Accepted core metadata changed in flight')
        c.need(c.historical_A1_ready(task) == old_reference, 'Original A1 reference metadata changed in flight')
        accepted = dict(status='extension_independently_verified', task_index=index, family=task['family'],
                        sample=task['sample'], budget=task['budget'], campaign_gate_sha256=gate_hash,
                        release_gate_sha256=gate['release_gate']['sha256'], worker_sha256=c.sha(__file__),
                        task_manifest_sha256=gate['task_manifest']['sha256'], dependencies=dependencies,
                        result=proof, artifact_stats=c.artifact_stats(proof['artifacts']),
                        job=attempt['job'], array_job=attempt['array_job'], array_task=index,
                        adapter_log_sha256=c.sha(unit / 'adapter.log'), completed_at=c.utc())
        c.write(unit / 'acceptance.json', accepted)
        (unit / 'EXTENSION_VERIFIED_COMPLETE').write_text(c.sha(unit / 'acceptance.json') + '\n')
        attempt.update(status='verified', completed_at=c.utc(), acceptance_sha256=c.sha(unit / 'acceptance.json'))
        c.write(unit / 'attempt.json', attempt)
        print(f'EXTENSION_VERIFIED_COMPLETE {index} {task["family"]}', flush=True)
    except BaseException as error:
        attempt.update(status='failed_preserved', error=repr(error), traceback=traceback.format_exc(),
                       failed_at=c.utc(), automatic_retry=False, artifacts_preserved=True)
        c.write(unit / 'attempt.json', attempt)
        raise


if __name__ == '__main__':
    main()
