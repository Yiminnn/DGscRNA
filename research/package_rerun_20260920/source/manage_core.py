"""Metadata-only, release-gated scheduler for the 726 packaged GBM core tasks.

No computation or evaluation is imported. Only the root starts this controller.
One array runs at a time; scientific failures require review, never a new method.
"""
from datetime import datetime, timedelta, timezone
from pathlib import Path
import argparse
import hashlib
import json
import os
import pwd
import re
import socket
import subprocess
import time
import uuid


HERE = Path(__file__).resolve().parent
ACTIVE = {'PENDING', 'RUNNING', 'CONFIGURING', 'COMPLETING', 'REQUEUED',
          'SUSPENDED', 'RESIZING', 'REQUEUE_FED', 'REQUEUE_HOLD'}
LIMITS = dict(queue_ceiling=900, queue_reserve=64, running_ceiling=256,
              array_concurrency=16, max_own_outstanding=64, max_active_arrays=1,
              poll_seconds=60)
CONFIG = dict(preset='gbm-reference', dataset='GSE274546', route='all',
              library='all', cutoff='all', seed=42, deg_workers=4, require_slurm=True)
REQUIRED_FILES = {'manage_core.py', 'core_task.py', 'core_array.sbatch',
                  'verify_packaged_parity.py', 'verify_packaged_r_artifacts.R',
                  'evaluate_core_lfine.py', 'evaluate_core_lfine_sources.json'}
PILOTS = {'anchor': ('TKU4163', 'hvg2000', 1),
          'full_roster': ('TKU4163', 'hvg2000', 192),
          'medium_all_genes': ('NL022', 'all', 1),
          'large_hvg2000': ('SN040', 'hvg2000', 1)}


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


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.part.' + uuid.uuid4().hex)
    with temporary.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def checked_file(record):
    path = Path(record['path']).resolve()
    need(path.is_file() and sha(path) == record['sha256'], f'Frozen artifact mismatch: {path}')
    return path


def gate_metadata(path):
    """Read only release/source/proof metadata; workers check numerical inputs."""
    gate = load(path)
    need(gate.get('status') == 'release_and_pilots_verified', 'Actual release gate is not ready')
    need(gate['config'] == CONFIG, 'Core scientific configuration changed')
    publication = gate['publication']
    need(publication['status'] == 'published' and publication['branch'] == 'align-r-reference',
         'Required authorized branch publication is incomplete')
    need(re.fullmatch('[0-9a-f]{40}', publication['commit']), 'Invalid publication commit')
    checked_file(publication['receipt'])
    checked_file(gate['wheel'])
    checked_file(gate['markers'])
    need(REQUIRED_FILES <= set(gate['launcher_files']), 'Gate must bind controller and worker/verifier launchers')
    for name, record in gate['launcher_files'].items():
        path_in_gate = checked_file(record)
        if name in REQUIRED_FILES:
            need(path_in_gate == HERE / name, f'Unexpected launcher location: {name}')
    tasks = load(checked_file(gate['task_manifest']))
    need(len(tasks) == 726 and [t['index'] for t in tasks] == list(range(726)), 'Incomplete task index roster')
    need(len({(t['sample'], t['budget']) for t in tasks}) == 726, 'Duplicate sample/budget tasks')
    root = Path(gate['output_root']).resolve()
    for task in tasks:
        need(re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*', task['sample']), 'Unsafe sample identifier')
        need(task['budget'] in {'hvg500', 'hvg1000', 'hvg2000', 'hvg3000', 'hvg5000', 'all'},
             'Unsupported task budget')
    pilots = {item['role']: item for item in gate['pilots']}
    need(len(pilots) == len(gate['pilots']) and set(pilots) == set(PILOTS), 'Missing/duplicate pilot roles')
    for role, (sample, budget, n_conditions) in PILOTS.items():
        proof_path = checked_file(pilots[role])
        proof = load(proof_path)
        need(proof.get('status') == 'passed_exact' and proof['sample'] == sample
             and proof['budget'] == budget and proof['terminal_conditions'] == n_conditions,
             f'Invalid pilot proof: {role}')
        plan_path = proof_path.parent / 'run_config.json'
        need(sha(plan_path) == proof['run_config_sha256'], f'Pilot configuration changed: {role}')
        need(load(plan_path)['runtime'] == gate['runtime']['fingerprint'], f'Pilot runtime differs: {role}')
    return gate, tasks, root


def command(argv, timeout=55):
    result = subprocess.run(argv, text=True, capture_output=True, timeout=timeout)
    need(result.returncode == 0, f'Metadata command failed: {argv}: {result.stderr[-2000:]}')
    return result.stdout


def queue():
    text = command(['squeue', '--array', '-u', 'yimin', '-h', '-o', '%i|%T|%100j'])
    rows = []
    for line in text.splitlines():
        fields = [x.strip() for x in line.split('|')]
        need(len(fields) == 3 and re.fullmatch(r'\d+(?:_\d+)?', fields[0]),
             f'Queue is not fully expanded: {line}')
        rows.append(dict(job=fields[0], state=fields[1], name=fields[2]))
    need(len({r['job'] for r in rows}) == len(rows), 'Duplicate queue identifiers')
    return rows


def accounting(jobs=None, *, name=None, since=None):
    if not jobs and not name:
        return {}
    argv = ['sacct', '--array', '-n', '-P', '-u', 'yimin',
            '--format=JobID%80,JobName%100,State%40,ExitCode,Elapsed,MaxRSS,ReqMem,NodeList,Reason']
    if jobs:
        argv += ['-j', ','.join(sorted(set(jobs)))]
    if name:
        # SLURM interprets naive timestamps in site time. Search from the prior
        # calendar day so a UTC/local offset cannot hide a just-submitted job.
        start = (datetime.fromisoformat(since) - timedelta(days=1)).strftime('%Y-%m-%d')
        argv += ['--name', name, '-S', start]
    rows = {}
    for line in command(argv).splitlines():
        fields = [x.strip() for x in line.split('|')]
        need(len(fields) >= 9, f'Malformed accounting row: {line}')
        job, job_name, state, exitcode, elapsed, rss, memory, nodes, reason = fields[:9]
        if not re.fullmatch(r'\d+(?:_\d+)?(?:\.[A-Za-z0-9_-]+)?', job):
            continue
        item = dict(job=job, name=job_name, state=state.split()[0].rstrip('+'),
                    exit_code=exitcode, elapsed=elapsed, max_rss=rss,
                    requested_memory=memory, nodes=nodes, reason=reason)
        need(job not in rows or rows[job] == item, f'Ambiguous accounting: {job}')
        rows[job] = item
    return rows


def output_for(root, task):
    return root / 'core' / task['sample'] / task['budget']


def acceptance(root, task, gate_hash, wheel_hash):
    output = output_for(root, task)
    flag = output / 'PACKAGE_VERIFIED_COMPLETE'
    if not flag.exists():
        return None
    path = output / 'package_acceptance.json'
    need(path.is_file() and flag.read_text().strip() == sha(path), f'Invalid acceptance flag: {output}')
    accepted = load(path)
    need(accepted.get('status') == 'package_run_independently_verified'
         and accepted['gate_sha256'] == gate_hash and accepted['wheel_sha256'] == wheel_hash,
         f'Acceptance belongs to a different release: {output}')
    need(all(accepted[k] == task[k] for k in ('sample', 'budget', 'index'))
         and accepted['terminal_conditions'] == 192, f'Incomplete/incorrect task acceptance: {output}')
    need(accepted['run_manifest_sha256'] == sha(output / 'run_manifest.json')
         and accepted['parity_sha256'] == sha(output / 'packaged_parity.json'),
         f'Accepted metadata changed: {output}')
    proof = load(output / 'packaged_parity.json')
    need(proof.get('status') == 'passed_exact' and proof['terminal_conditions'] == 192
         and proof['sample'] == task['sample'] and proof['budget'] == task['budget'],
         f'Independent parity did not pass: {output}')
    need(proof['run_manifest_sha256'] == accepted['run_manifest_sha256']
         and proof['run_config_sha256'] == sha(output / 'run_config.json'),
         f'Parity is not bound to this run: {output}')
    evaluation_path = Path(accepted['lfine_manifest']).resolve()
    need(evaluation_path == (root / 'evaluation/units' / task['sample'] / task['budget'] / 'manifest.json').resolve(),
         'Lfine manifest differs from the canonical unit evaluation path')
    need(accepted['lfine_threshold_rows'] == accepted['lfine_valid_threshold_rows'] == 384,
         f'Incomplete Lfine threshold rows: {output}')
    need(sha(evaluation_path) == accepted['lfine_manifest_sha256']
         and (evaluation_path.parent / 'COMPLETE').read_text().strip() == accepted['lfine_manifest_sha256'],
         f'Lfine manifest/flag changed: {output}')
    evaluation = load(evaluation_path)
    metrics_hash = sha(evaluation_path.parent / 'metrics.csv.gz')
    need(evaluation.get('status') == 'completed' and evaluation['n_conditions'] == evaluation['n_valid'] == 384
         and evaluation['sample'] == task['sample'] and evaluation['budget'] == task['budget']
         and evaluation['task_index'] == task['index'], f'Lfine evaluation is incomplete or mismatched: {output}')
    need(metrics_hash == accepted['lfine_metrics_sha256'] == evaluation['outputs']['metrics.csv.gz'],
         f'Lfine metric payload changed: {output}')
    need(evaluation['script_sha256'] == sha(HERE / 'evaluate_core_lfine.py')
         and evaluation['frozen_semantics'] == load(HERE / 'evaluate_core_lfine_sources.json'),
         f'Lfine evaluation used different frozen source/semantics: {output}')
    return dict(path=str(path), sha256=sha(path), parity_sha256=accepted['parity_sha256'],
                run_manifest_sha256=accepted['run_manifest_sha256'], index=task['index'],
                sample=task['sample'], budget=task['budget'], job=accepted['job'],
                array_job=accepted['array_job'], array_task=accepted['array_task'],
                lfine_manifest=str(evaluation_path), lfine_manifest_sha256=accepted['lfine_manifest_sha256'],
                lfine_metrics_sha256=metrics_hash, lfine_threshold_rows=384)


def registry(state, gate_hash, batch_size):
    need(state['schema'] == 1 and state['gate_sha256'] == gate_hash
         and state['controller_sha256'] == sha(__file__) and state['limits'] == LIMITS
         and state['batch_size'] == batch_size, 'Controller state/gate/options changed')
    assigned, jobs = set(), set()
    for submission in state['submissions']:
        if submission['status'] != 'submitted':
            continue
        need(re.fullmatch(r'\d+', submission['job_id']) and submission['job_id'] not in jobs,
             'Invalid/duplicate registered array job')
        jobs.add(submission['job_id'])
        indices = submission['indices']
        need(0 < len(indices) <= batch_size and len(set(indices)) == len(indices)
             and all(isinstance(i, int) and 0 <= i < 726 for i in indices), 'Invalid array indices')
        need(not assigned.intersection(indices), 'Duplicate task submission requires explicit review')
        need(submission['concurrency'] == LIMITS['array_concurrency'], 'Changed array concurrency')
        assigned.update(indices)
    return assigned, jobs


def matching_intent(intent, queued):
    matches = {r['job'].split('_')[0] for r in queued if r['name'] == intent['name']}
    records = accounting(name=intent['name'], since=intent['created_at'])
    matches.update(r['job'].split('_')[0].split('.')[0] for r in records.values()
                   if r['name'] == intent['name'])
    need(len(matches) <= 1, 'Submission name matched multiple jobs; inspect before proceeding')
    return next(iter(matches), None)


def recover_intent(state, control, queued):
    intent = state.get('submission_intent')
    if not intent:
        return True
    job = matching_intent(intent, queued)
    if job:
        state['submissions'].append(dict(**intent, status='submitted', job_id=job,
                                         recovered_at=utc(), submission_reply_reconciled=True))
        state.pop('submission_intent')
        write(control / 'state.json', state)
        return True
    # Only an explicit scheduler submit-limit rejection proves no worker was
    # created. Timeouts/transport errors remain unresolved rather than duplicate.
    if intent.get('explicit_capacity_rejection'):
        state['submissions'].append(dict(**intent, status='rejected_not_submitted'))
        state['retry_after_epoch'] = time.time() + 120
        state.pop('submission_intent')
        write(control / 'state.json', state)
        return True
    return False


def publish(control, state, status, **fields):
    record = dict(status=status, updated_at=utc(), gate_sha256=state['gate_sha256'],
                  controller_sha256=state['controller_sha256'], limits=LIMITS,
                  verified_tasks=len(state['accepted']), expected_tasks=726,
                  submitted_arrays=[s['job_id'] for s in state['submissions'] if s['status'] == 'submitted'],
                  scientific_computation_performed_by_controller=False,
                  automatic_scientific_retry=False, automatic_PTC_launch=False, **fields)
    write(control / 'status.json', record)
    print(json.dumps({k: record[k] for k in ['status', 'updated_at', 'verified_tasks', 'expected_tasks']}), flush=True)


def task_failures(state, tasks, root, records, queued, accepted):
    current = {row['job']: row for row in queued}
    busy, failures = False, []
    for submission in state['submissions']:
        if submission['status'] != 'submitted':
            continue
        for index in submission['indices']:
            job = f"{submission['job_id']}_{index}"
            record = records.get(job)
            if job in current or record is None or record['state'] in ACTIVE:
                busy = True
                continue
            failed = record['state'] != 'COMPLETED' or record['exit_code'] != '0:0'
            if failed or str(index) not in accepted:
                output = output_for(root, tasks[index])
                audits = sorted((root / 'control/attempts').glob(f'{index:04d}_*.json'))
                failures.append(dict(index=index, sample=tasks[index]['sample'], budget=tasks[index]['budget'],
                                     scheduler=record,
                                     step_accounting=[r for name, r in records.items() if name.startswith(job + '.')],
                                     output=str(output), output_exists=output.exists(),
                                     worker_attempt_receipts=[str(p) for p in audits],
                                     reason='scheduler_failure' if failed else 'missing_verified_acceptance',
                                     retry_policy='No automatic retry; preserve artifacts and review the same frozen protocol'))
    return busy, failures


def may_submit(queued, n):
    running = sum(row['state'] not in {'PENDING', 'REQUEUED', 'REQUEUE_HOLD'} for row in queued)
    return (len(queued) + n + LIMITS['queue_reserve'] < LIMITS['queue_ceiling']
            and running + LIMITS['array_concurrency'] < LIMITS['running_ceiling'])


def submit(state, control, gate_path, gate, indices):
    # Read the queue again immediately before the atomic intent and sbatch call.
    queued = queue()
    if not may_submit(queued, len(indices)):
        return False
    name = 'pkgcore_' + state['gate_sha256'][:10] + '_' + uuid.uuid4().hex[:10]
    array = ','.join(map(str, indices)) + '%' + str(LIMITS['array_concurrency'])
    argv = ['sbatch', '--parsable', '--job-name=' + name, '--array=' + array,
            gate['launcher_files']['core_array.sbatch']['path'], str(gate_path)]
    intent = dict(id=uuid.uuid4().hex, name=name, created_at=utc(), indices=indices,
                  concurrency=LIMITS['array_concurrency'], argv=argv,
                  queued_before=len(queued), gate_sha256=state['gate_sha256'])
    state['submission_intent'] = intent
    write(control / 'state.json', state)
    try:
        result = subprocess.run(argv, text=True, capture_output=True, timeout=45)
        intent.update(returncode=result.returncode, stdout=result.stdout, stderr=result.stderr,
                      reply_at=utc())
        response = result.stdout.strip()
        if result.returncode == 0 and re.fullmatch(r'\d+(?:;[^\s;]+)?', response):
            state['submissions'].append(dict(**intent, status='submitted', job_id=response.split(';')[0]))
            state.pop('submission_intent')
        else:
            intent['explicit_capacity_rejection'] = bool(
                result.returncode != 0 and not response and
                re.search(r'Batch job submission failed:.*(?:AssocMaxSubmitJobLimit|'
                          r'QOSMaxSubmitJobPerUserLimit|MaxSubmitJobs|[Jj]ob submit limit)', result.stderr))
    except (OSError, subprocess.TimeoutExpired) as error:
        intent.update(submission_exception=repr(error), reply_at=utc(),
                      explicit_capacity_rejection=False)
    write(control / 'state.json', state)
    # Preserve scheduler metadata even when another controller consumed headroom.
    after = queue()
    state['last_post_submission_queue'] = dict(at=utc(), expanded=len(after),
                                             running=sum(r['state'] == 'RUNNING' for r in after))
    write(control / 'state.json', state)
    need(len(after) < LIMITS['queue_ceiling'] and sum(r['state'] == 'RUNNING' for r in after) < 256,
         'Concurrent scheduler activity exhausted limits; no jobs cancelled, further submissions stopped')
    return True


def iteration(state, control, gate_path, gate_hash, batch_size):
    need(sha(gate_path) == gate_hash, 'Release gate changed while controller was running')
    gate, tasks, root = gate_metadata(gate_path)
    assigned, job_ids = registry(state, gate_hash, batch_size)
    queued = queue()
    if not recover_intent(state, control, queued):
        publish(control, state, 'submission_reply_uncertain', submission_intent=state['submission_intent'])
        return False
    assigned, job_ids = registry(state, gate_hash, batch_size)
    prefix = 'pkgcore_' + gate_hash[:10] + '_'
    unknown = [r for r in queued if r['name'].startswith(prefix) and r['job'].split('_')[0] not in job_ids]
    need(not unknown, f'Unregistered jobs claim this gate; adoption review required: {unknown}')
    own_running = [r for r in queued if r['job'].split('_')[0] in job_ids
                   and r['state'] not in {'PENDING', 'REQUEUED', 'REQUEUE_HOLD'}]
    need(len(own_running) <= LIMITS['array_concurrency'],
         'Registered campaign running-task cap exceeded; stop submissions without cancelling jobs')
    records = accounting(job_ids)
    for task in tasks:
        proof = acceptance(root, task, gate_hash, gate['wheel']['sha256'])
        key = str(task['index'])
        if proof is None:
            need(key not in state['accepted'], f'Previously accepted receipt disappeared: {key}')
            continue
        need(task['index'] in assigned, f'Accepted task has no registered submission: {key}')
        submission = next(s for s in state['submissions']
                          if s['status'] == 'submitted' and task['index'] in s['indices'])
        need(str(proof['array_job']) == submission['job_id'] and int(proof['array_task']) == task['index'],
             f'Acceptance job differs from registered worker: {key}')
        need(key not in state['accepted'] or state['accepted'][key] == proof,
             f'Previously accepted receipt changed: {key}')
        state['accepted'][key] = proof
    write(control / 'state.json', state)
    busy, failures = task_failures(state, tasks, root, records, queued, state['accepted'])
    if failures:
        write(control / 'failures.json', dict(at=utc(), failures=failures,
                                             no_jobs_cancelled=True, automatic_retry=False))
        publish(control, state, 'blocked_failed_tasks', failures=failures)
        return True
    outstanding = len(assigned - {int(i) for i in state['accepted']})
    need(outstanding <= LIMITS['max_own_outstanding'], 'Own outstanding task cap exceeded')
    if len(state['accepted']) == 726 and not busy:
        receipt = dict(status='all_726_packaged_GBM_core_tasks_independently_verified',
                       completed_at=utc(), gate_sha256=gate_hash, controller_sha256=sha(__file__),
                       task_manifest_sha256=gate['task_manifest']['sha256'],
                       accepted_tasks=726, terminal_conditions=726 * 192,
                       lfine_valid_threshold_rows=726 * 384,
                       acceptances=state['accepted'], PTC_started=False)
        target = root / 'GBM_CORE_COMPLETE'
        if target.exists():
            prior = load(target)
            need(prior['gate_sha256'] == gate_hash and prior['acceptances'] == state['accepted'],
                 'Existing completion receipt differs')
        else:
            write(target, receipt)
        publish(control, state, 'completed', completion_receipt=str(target))
        return True
    if busy:
        publish(control, state, 'waiting_active_array', own_outstanding=outstanding,
                expanded_queue=len(queued), running=sum(r['state'] == 'RUNNING' for r in queued))
        return False
    if time.time() < state.get('retry_after_epoch', 0):
        publish(control, state, 'waiting_scheduler_submission_cooldown')
        return False
    ready = [t['index'] for t in tasks if t['index'] not in assigned]
    need(ready, 'Incomplete campaign has no schedulable task; inspect receipts')
    indices = ready[:min(batch_size, LIMITS['max_own_outstanding'] - outstanding)]
    # Never launch a scientific retry into a preserved partial output directory.
    unexpected = [str(output_for(root, tasks[i])) for i in indices if output_for(root, tasks[i]).exists()]
    need(not unexpected, f'Unregistered preserved outputs require review: {unexpected}')
    sent = bool(indices) and may_submit(queued, len(indices)) and submit(state, control, gate_path, gate, indices)
    publish(control, state, 'submitted' if sent else 'waiting_queue_capacity',
            own_outstanding=outstanding + (len(indices) if sent else 0),
            expanded_queue=len(queued), queue_reserve=LIMITS['queue_reserve'])
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gate', required=True, type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--once', action='store_true')
    mode.add_argument('--watch', action='store_true')
    parser.add_argument('--batch-size', type=int, choices=[32, 64], default=64)
    args = parser.parse_args()
    need(pwd.getpwuid(os.getuid()).pw_name == 'yimin', 'Only the authorized local yimin identity may run this controller')
    need(os.environ.get('SLURM_JOB_ID'), 'Start this administrative controller within the existing SLURM allocation')
    need(not os.sys.flags.optimize, 'Do not start the controller with Python -O')
    gate_path = args.gate.resolve()
    gate, tasks, root = gate_metadata(gate_path)
    gate_hash = sha(gate_path)
    control = root / 'control/core_manager'
    control.mkdir(parents=True, exist_ok=True)
    lock = control / 'manager.lock'
    try:
        lock.mkdir()
    except FileExistsError as error:
        raise RuntimeError(f'Controller lock exists: {lock}; inspect owner/job before restart, never clear a live lock') from error
    state = None
    try:
        write(lock / 'owner.json', dict(job=os.environ['SLURM_JOB_ID'],
                                       host=socket.gethostname(), pid=os.getpid(), started_at=utc(),
                                       gate_sha256=gate_hash, controller_sha256=sha(__file__)))
        path = control / 'state.json'
        if path.exists():
            state = load(path)
            registry(state, gate_hash, args.batch_size)
        else:
            state = dict(schema=1, created_at=utc(), gate_sha256=gate_hash,
                         controller_sha256=sha(__file__), limits=LIMITS, batch_size=args.batch_size,
                         submissions=[], accepted={})
            write(path, state)
        while True:
            if (control / 'STOP').exists():
                publish(control, state, 'stopped_by_local_flag', jobs_left_untouched=True)
                break
            finished = iteration(state, control, gate_path, gate_hash, args.batch_size)
            if finished or args.once:
                break
            time.sleep(LIMITS['poll_seconds'])
    except BaseException as error:
        if state is not None:
            publish(control, state, 'blocked_controller_review', error=repr(error), jobs_left_untouched=True)
        raise
    finally:
        (lock / 'owner.json').unlink(missing_ok=True)
        lock.rmdir()


if __name__ == '__main__':
    main()
