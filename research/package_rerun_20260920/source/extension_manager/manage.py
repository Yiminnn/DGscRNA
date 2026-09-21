"""Durable metadata-only scheduler for dependency-ready packaged GBM extensions."""
from pathlib import Path
import argparse
import json
import os
import re
import socket
import subprocess
import time
import uuid
import common as c
import receipts

queue, accounting = c.core.queue, c.core.accounting
ACTIVE = c.core.ACTIVE


def registry(state, gate_hash, tasks, batch_size):
    c.need(state['schema'] == 1 and state['campaign_gate_sha256'] == gate_hash
           and state['controller_sha256'] == c.sha(__file__) and state['limits'] == c.LIMITS
           and state['batch_size'] == batch_size and state['task_count'] == len(tasks), 'Controller state/options/source changed')
    assigned, jobs = set(), set()
    for submission in state['submissions']:
        if submission['status'] != 'submitted':
            continue
        job, indices = submission['job_id'], submission['indices']
        c.need(re.fullmatch(r'\d+', job) and job not in jobs, 'Duplicate/invalid registered array')
        c.need(0 < len(indices) <= batch_size and len(set(indices)) == len(indices)
               and all(type(i) is int and 0 <= i < len(tasks) for i in indices), 'Invalid extension indices')
        c.need(not assigned.intersection(indices) and submission['concurrency'] == c.LIMITS['initial_concurrency'],
               'Duplicate scientific submission or changed concurrency')
        assigned.update(indices); jobs.add(job)
    return assigned, jobs


def may_submit(queued, n):
    running = sum(r['state'] not in {'PENDING', 'REQUEUED', 'REQUEUE_HOLD'} for r in queued)
    return len(queued) + n + c.LIMITS['queue_reserve'] < c.LIMITS['queue_ceiling'] and \
        running + c.LIMITS['initial_concurrency'] + c.LIMITS['core_running_reserve'] < c.LIMITS['running_ceiling']


def throttle_target(queued, own_job):
    running = [r for r in queued if r['state'] not in {'PENDING', 'REQUEUED', 'REQUEUE_HOLD'}]
    own_running = sum(r['job'].split('_')[0] == own_job for r in running)
    unrelated = len(running) - own_running
    available = c.LIMITS['running_ceiling'] - 1 - c.LIMITS['core_running_reserve'] - unrelated
    target = max(c.LIMITS['minimum_concurrency'], min(c.LIMITS['maximum_concurrency'], available))
    return dict(target=target, available_with_core_reserve=available, all_user_running=len(running),
                own_running=own_running, unrelated_running=unrelated,
                reserve_temporarily_unavailable=available < c.LIMITS['minimum_concurrency'])


def adapt_throttle(state, control, queued):
    active = [s for s in state['submissions'] if s['status'] == 'submitted'
              and any(r['job'].split('_')[0] == s['job_id'] for r in queued)]
    c.need(len(active) <= 1, 'More than one active extension array')
    if not active:
        return
    submission = active[0]
    # Recheck global state immediately before a scheduling-only change.
    fresh = queue()
    decision = throttle_target(fresh, submission['job_id'])
    current = submission.get('current_throttle', submission['concurrency'])
    decision.update(at=c.utc(), array_job=submission['job_id'], previous=current,
                    scheduling_only=True, scientific_parameters_unchanged=True,
                    unrelated_controllers_not_constrained=True)
    pending = state.get('throttle_intent')
    if decision['target'] != current or pending:
        # A prior scontrol may have succeeded before its reply/state update.
        # Never infer the live throttle from cached state after such a crash.
        # Setting the newly computed target again is safe and idempotent.
        if pending:
            decision['reconciles_unresolved_intent'] = pending
        argv = ['scontrol', 'update', 'JobId=' + submission['job_id'],
                'ArrayTaskThrottle=' + str(decision['target'])]
        decision['argv'] = argv
        state['throttle_intent'] = decision; c.write(control / 'state.json', state)
        result = subprocess.run(argv, text=True, capture_output=True, timeout=45)
        decision.update(returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
        state.setdefault('throttle_decisions', []).append(decision)
        c.write(control / 'state.json', state)
        c.need(result.returncode == 0, 'Own-array throttle update failed; preserve jobs and review')
        submission['current_throttle'] = decision['target']
        state.pop('throttle_intent', None)
    state['last_capacity_decision'] = decision
    c.write(control / 'state.json', state)


def recover_intent(state, control, queued):
    intent = state.get('submission_intent')
    if not intent:
        return True
    matching = {r['job'].split('_')[0] for r in queued if r['name'] == intent['name']}
    history = accounting(name=intent['name'], since=intent['created_at'])
    matching.update(r['job'].split('_')[0].split('.')[0] for r in history.values() if r['name'] == intent['name'])
    c.need(len(matching) <= 1, 'Submission intent matched multiple jobs; review required')
    if matching:
        state['submissions'].append(dict(**intent, status='submitted', job_id=next(iter(matching)),
                                         recovered_at=c.utc(), reply_reconciled=True))
        state.pop('submission_intent'); c.write(control / 'state.json', state)
        return True
    if intent.get('explicit_capacity_rejection'):
        state['submissions'].append(dict(**intent, status='rejected_not_submitted'))
        state.pop('submission_intent'); state['retry_after_epoch'] = time.time() + 120
        c.write(control / 'state.json', state)
        return True
    return False


def publish(control, state, status, **fields):
    value = dict(status=status, updated_at=c.utc(), campaign_gate_sha256=state['campaign_gate_sha256'],
                 controller_sha256=state['controller_sha256'], limits=c.LIMITS,
                 verified_tasks=len(state['accepted']), expected_tasks=state['task_count'],
                 automatic_scientific_retry=False, scientific_computation_by_controller=False,
                 unrelated_jobs_changed=False, PTC_started=False,
                 capacity=state.get('last_capacity_decision'),
                 submitted_arrays=[s['job_id'] for s in state['submissions'] if s['status'] == 'submitted'], **fields)
    c.write(control / 'status.json', value)
    print(json.dumps({k:value[k] for k in ['status', 'updated_at', 'verified_tasks', 'expected_tasks']}), flush=True)


artifact_stats = c.artifact_stats


def acceptance(root, task, gate_hash, release_hash, previous=None, *, final=False):
    unit = c.worker_receipt_root(root, task['index'])
    if not (unit / 'EXTENSION_VERIFIED_COMPLETE').exists():
        return None
    value = c.receipt(unit, 'acceptance.json', 'EXTENSION_VERIFIED_COMPLETE')
    c.need(value['status'] == 'extension_independently_verified'
           and value['campaign_gate_sha256'] == gate_hash and value['release_gate_sha256'] == release_hash
           and value['worker_sha256'] == c.sha(c.HERE / 'worker.py')
           and value['task_index'] == task['index'] and value['family'] == task['family']
           and value['sample'] == task['sample'] and value['budget'] == task['budget'],
           'Worker acceptance belongs to another task/release')
    if previous is not None and not final:
        c.need(previous['sha256'] == c.sha(unit / 'acceptance.json'), 'Worker acceptance changed')
        c.need(previous['artifact_stats'] == artifact_stats(value['result']['artifacts']),
               'An accepted artifact stat changed; stop for integrity review')
        return previous
    c.need(value['artifact_stats'] == artifact_stats(value['result']['artifacts']),
           'An artifact changed after the worker deep-checksum audit')
    proof = receipts.validate(task, release_hash, deep=final)
    c.need(value['result'] == proof, 'Worker accepted artifacts changed')
    return dict(path=str(unit / 'acceptance.json'), sha256=c.sha(unit / 'acceptance.json'),
                index=task['index'], family=task['family'], sample=task['sample'], budget=task['budget'],
                array_job=str(value['array_job']), array_task=int(value['array_task']), job=str(value['job']),
                terminal_conditions=proof['terminal_conditions'], lfine_rows=proof['lfine_rows'],
                artifact_stats=artifact_stats(proof['artifacts']))


def failures_and_busy(state, tasks, root, records, queued):
    current = {r['job'] for r in queued}
    busy, failures = False, []
    for submission in state['submissions']:
        if submission['status'] != 'submitted':
            continue
        for index in submission['indices']:
            job = f'{submission["job_id"]}_{index}'
            row = records.get(job)
            if job in current or row is None or row['state'] in ACTIVE:
                busy = True
                continue
            if row['state'] != 'COMPLETED' or row['exit_code'] != '0:0' or str(index) not in state['accepted']:
                attempt = c.worker_receipt_root(root, index) / 'attempt.json'
                failures.append(dict(index=index, family=tasks[index]['family'], scheduler=row,
                                     steps=[v for k,v in records.items() if k.startswith(job + '.')],
                                     worker_attempt=str(attempt), worker_metadata=c.load(attempt) if attempt.exists() else None,
                                     output=tasks[index]['output'], reason='scheduler_failure_or_missing_verified_receipt',
                                     automatic_retry=False, artifacts_preserved=True))
    return busy, failures


def submit(state, control, gate_path, gate, tasks, release, root, indices):
    # Repeat source/dependency checks and queue admission immediately before sbatch.
    c.gate_metadata(gate_path)
    cache = {}
    c.need(all(c.dependency_proofs(tasks[i], release, root, cache) is not None for i in indices),
           'Core dependency disappeared before submission')
    c.need(all(not Path(tasks[i]['output']).exists() and not c.worker_receipt_root(root, i).exists() for i in indices),
           'Preserved output/worker attempt exists; do not retry automatically')
    c.need(all(c.historical_A1_ready(tasks[i]) is not None for i in indices), 'Original A1 reference is incomplete')
    c.need(all(all(str(dep) in state['accepted'] for dep in tasks[i]['extension_dependencies']) for i in indices),
           'Required noDR export/verification has not completed')
    profiles = {tasks[i]['profile'] for i in indices}
    c.need(len(profiles) == 1, 'One array must use one resource profile')
    profile_name = next(iter(profiles)); profile = c.PROFILES[profile_name]
    queued = queue()
    if not may_submit(queued, len(indices)):
        return False
    name = 'pkgext_' + state['campaign_gate_sha256'][:10] + '_' + uuid.uuid4().hex[:10]
    argv = ['sbatch', '--parsable', '--job-name=' + name,
            '--array=' + ','.join(map(str, indices)) + '%' + str(c.LIMITS['initial_concurrency']),
            '--cpus-per-task=' + str(profile['cpus']), '--mem=' + profile['memory'], '--time=' + profile['time'],
            gate['launchers']['array.sbatch']['path'], str(gate_path)]
    intent = dict(id=uuid.uuid4().hex, name=name, created_at=c.utc(), indices=indices,
                  concurrency=c.LIMITS['initial_concurrency'], argv=argv, queued_before=len(queued),
                  campaign_gate_sha256=state['campaign_gate_sha256'], resource_profile=profile_name)
    state['submission_intent'] = intent; c.write(control / 'state.json', state)
    try:
        result = subprocess.run(argv, text=True, capture_output=True, timeout=45)
        response = result.stdout.strip()
        intent.update(returncode=result.returncode, stdout=result.stdout, stderr=result.stderr, reply_at=c.utc())
        if result.returncode == 0 and re.fullmatch(r'\d+(?:;[^\s;]+)?', response):
            state['submissions'].append(dict(**intent, status='submitted', job_id=response.split(';')[0]))
            state.pop('submission_intent')
        else:
            intent['explicit_capacity_rejection'] = bool(result.returncode != 0 and not response and re.search(
                r'Batch job submission failed:.*(?:AssocMaxSubmitJobLimit|QOSMaxSubmitJobPerUserLimit|MaxSubmitJobs|[Jj]ob submit limit)',
                result.stderr))
    except (OSError, subprocess.TimeoutExpired) as error:
        intent.update(submission_exception=repr(error), explicit_capacity_rejection=False, reply_at=c.utc())
    c.write(control / 'state.json', state)
    after = queue()
    state['last_post_submission_queue'] = dict(at=c.utc(), expanded=len(after), running=sum(r['state'] == 'RUNNING' for r in after))
    c.write(control / 'state.json', state)
    c.need(len(after) < c.LIMITS['queue_ceiling'] and sum(r['state'] == 'RUNNING' for r in after) < 256,
           'Concurrent campaigns exhausted queue/run headroom; preserve all jobs and stop submissions')
    return True


def iteration(state, control, gate_path, gate_hash, batch_size):
    c.need(c.sha(gate_path) == gate_hash, 'Campaign gate changed during execution')
    gate, tasks, release, root = c.gate_metadata(gate_path)
    assigned, jobs = registry(state, gate_hash, tasks, batch_size)
    queued = queue()
    if not recover_intent(state, control, queued):
        publish(control, state, 'submission_reply_uncertain', submission_intent=state['submission_intent'])
        return False
    assigned, jobs = registry(state, gate_hash, tasks, batch_size)
    prefix = 'pkgext_' + gate_hash[:10] + '_'
    c.need(not [r for r in queued if r['name'].startswith(prefix) and r['job'].split('_')[0] not in jobs],
           'Unregistered extension jobs require explicit adoption review')
    own = [r for r in queued if r['job'].split('_')[0] in jobs]
    own_running = [r for r in own if r['state'] not in {'PENDING', 'REQUEUED', 'REQUEUE_HOLD'}]
    c.need(len(own) <= 64 and len(own_running) <= 64, 'Own outstanding/concurrency cap exceeded')
    records = accounting(jobs)
    for index in assigned:
        task = tasks[index]
        proof = acceptance(root, task, gate_hash, gate['release_gate']['sha256'], state['accepted'].get(str(index)))
        if proof is None:
            c.need(str(index) not in state['accepted'], 'A previously accepted worker receipt disappeared')
            continue
        sub = next(s for s in state['submissions'] if s['status'] == 'submitted' and index in s['indices'])
        c.need(proof['array_job'] == sub['job_id'] and proof['array_task'] == index, 'Acceptance worker differs from registered job')
        c.need(str(index) not in state['accepted'] or state['accepted'][str(index)] == proof, 'Previously accepted result changed')
        state['accepted'][str(index)] = proof
    c.write(control / 'state.json', state)
    busy, failures = failures_and_busy(state, tasks, root, records, queued)
    if failures:
        c.write(control / 'failures.json', dict(at=c.utc(), failures=failures, automatic_retry=False, jobs_left_untouched=True))
        publish(control, state, 'blocked_failed_tasks', failures=failures)
        return True
    outstanding = len(assigned - {int(k) for k in state['accepted']})
    c.need(outstanding <= 64, 'Outstanding cap exceeded')
    if len(state['accepted']) == len(tasks) and not busy:
        publish(control, state, 'final_integrity_audit', audited_tasks=0)
        for i, task in enumerate(tasks):
            proof = acceptance(root, task, gate_hash, gate['release_gate']['sha256'], final=True)
            c.need(proof == state['accepted'][str(task['index'])], 'Final accepted artifact audit differs')
            if (i + 1) % 16 == 0:
                publish(control, state, 'final_integrity_audit', audited_tasks=i + 1)
        completed = dict(status='all_declared_GBM_extensions_verified', campaign_gate_sha256=gate_hash,
                         completed_at=c.utc(), accepted_tasks=len(tasks), families=gate['families'],
                         terminal_conditions=sum(p['terminal_conditions'] for p in state['accepted'].values()),
                         lfine_rows=sum(p['lfine_rows'] for p in state['accepted'].values()),
                         acceptances=state['accepted'], optional_families=gate['optional_families'], PTC_started=False)
        target = root / ('GBM_EXTENSIONS_COMPLETE' if gate['scope'] == 'full_2617' else 'GBM_INITIAL_EXTENSIONS_COMPLETE')
        if target.exists():
            prior = c.load(target)
            c.need(prior['campaign_gate_sha256'] == gate_hash and prior['acceptances'] == state['accepted'],
                   'Existing completion receipt differs')
        else:
            c.write(target, completed)
        publish(control, state, 'completed', completion_receipt=str(target))
        return True
    if busy:
        adapt_throttle(state, control, queued)
        publish(control, state, 'waiting_active_array', own_outstanding=outstanding, expanded_queue=len(queued))
        return False
    if time.time() < state.get('retry_after_epoch', 0):
        publish(control, state, 'waiting_submission_cooldown')
        return False
    ready, waiting, cache = [], [], {}
    for task in tasks:
        if task['index'] in assigned:
            continue
        if c.dependency_proofs(task, release, root, cache) is None or \
           any(str(dep) not in state['accepted'] for dep in task['extension_dependencies']) or \
           c.historical_A1_ready(task) is None:
            waiting.append(task['index'])
        else:
            ready.append(task['index'])
    if not ready:
        publish(control, state, 'waiting_core_dependencies', waiting_tasks=len(waiting), ready_tasks=0)
        return False
    profile = tasks[ready[0]]['profile']
    indices = [i for i in ready if tasks[i]['profile'] == profile][:min(batch_size, 64 - outstanding)]
    sent = bool(indices) and may_submit(queued, len(indices)) and submit(state, control, gate_path, gate, tasks, release, root, indices)
    publish(control, state, 'submitted' if sent else 'waiting_queue_capacity', ready_tasks=len(ready),
            waiting_core_tasks=len(waiting), own_outstanding=outstanding + (len(indices) if sent else 0), expanded_queue=len(queued))
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', required=True, type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--once', action='store_true'); mode.add_argument('--watch', action='store_true')
    parser.add_argument('--batch-size', type=int, choices=[32, 64], default=64)
    args = parser.parse_args(); c.safe_identity()
    path = args.campaign.resolve(); gate_hash = c.sha(path)
    gate, tasks, release, root = c.gate_metadata(path)
    control = root / 'control/extension_manager'; control.mkdir(parents=True, exist_ok=True)
    lock = control / 'manager.lock'
    try:
        lock.mkdir()
    except FileExistsError as error:
        raise RuntimeError(f'Controller lock exists: {lock}; inspect owner/job before any restart') from error
    state = None
    try:
        c.write(lock / 'owner.json', dict(job=os.environ['SLURM_JOB_ID'], host=socket.gethostname(), pid=os.getpid(),
                                         started_at=c.utc(), campaign_gate_sha256=gate_hash, controller_sha256=c.sha(__file__)))
        if (control / 'state.json').exists():
            state = c.load(control / 'state.json'); registry(state, gate_hash, tasks, args.batch_size)
        else:
            state = dict(schema=1, created_at=c.utc(), campaign_gate_sha256=gate_hash, controller_sha256=c.sha(__file__),
                         limits=c.LIMITS, task_count=len(tasks), batch_size=args.batch_size, submissions=[], accepted={})
            c.write(control / 'state.json', state)
        while True:
            if (control / 'STOP').exists():
                publish(control, state, 'stopped_by_local_flag', jobs_left_untouched=True); break
            if iteration(state, control, path, gate_hash, args.batch_size) or args.once:
                break
            time.sleep(c.LIMITS['poll_seconds'])
    except BaseException as error:
        if state is not None:
            publish(control, state, 'blocked_controller_review', error=repr(error), jobs_left_untouched=True)
        raise
    finally:
        (lock / 'owner.json').unlink(missing_ok=True); lock.rmdir()


if __name__ == '__main__':
    main()
