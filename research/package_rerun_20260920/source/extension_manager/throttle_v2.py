"""Read-back verified SLURM array throttles; finished members are not fit failures."""
import re
import subprocess


WAITING = {'PENDING', 'REQUEUED', 'REQUEUE_HOLD'}


def own_rows(rows, job):
    return [r for r in rows if r['job'].split('_')[0] == job]


def finished_members_only(stderr, job):
    lines = [line.strip() for line in stderr.splitlines() if line.strip()]
    pattern = re.compile(re.escape(job) + r'_[0-9,\[\]%-]+: Job has already finished')
    return bool(lines) and all(pattern.fullmatch(line) for line in lines)


def readback(submission):
    result = subprocess.run(['scontrol', 'show', 'job', submission['job_id'], '-o'],
        text=True, capture_output=True, timeout=45)
    if result.returncode:
        raise RuntimeError('Cannot verify own array throttle: ' + result.stderr)
    records = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        fields = dict(word.split('=', 1) for word in line.split() if '=' in word)
        if fields.get('ArrayJobId') != submission['job_id'] or \
           fields.get('JobName') != submission['name'] or \
           not fields.get('UserId', '').startswith('yimin('):
            raise RuntimeError('Throttle readback includes an unowned array')
        value = fields.get('ArrayTaskThrottle', '')
        if not value.isdigit():
            raise RuntimeError('Throttle readback is missing its numeric cap')
        records.append(dict(job=fields['JobId'], array_task=fields.get('ArrayTaskId'),
            state=fields['JobState'], throttle=int(value)))
    if not records:
        raise RuntimeError('Empty own-array throttle readback')
    return records


def retire_completed_intent(state, control, records, queued, common):
    """Resolve an obsolete intent only using completed/accepted array evidence."""
    c = common
    pending = state.get('throttle_intent')
    if not pending:
        return
    job = pending['array_job']
    matches = [s for s in state['submissions'] if s.get('job_id') == job and s['status'] == 'submitted']
    c.need(len(matches) == 1 and not own_rows(queued, job), 'Unresolved throttle array is not terminal')
    if pending.get('returncode') not in (None, 0):
        c.need(finished_members_only(pending.get('stderr', ''), job), 'Unknown terminal-array throttle error requires review')
    rows = []
    for index in matches[0]['indices']:
        record = records.get(f'{job}_{index}')
        c.need(record is not None and record['state'] == 'COMPLETED' and record['exit_code'] == '0:0'
               and str(index) in state['accepted'], 'Cannot retire throttle before successful scientific acceptance/accounting')
        rows.append(record)
    decision = dict(at=c.utc(), array_job=job, action='superseded_terminal_array',
        superseded_unresolved_intent=pending, successful_throttle_update_claimed=False,
        completed_accounting=rows, scheduling_only=True, scientific_parameters_unchanged=True)
    state.setdefault('throttle_decisions', []).append(decision)
    state['last_capacity_decision'] = decision
    state.pop('throttle_intent')
    c.write(control / 'state.json', state)


def reconcile(state, control, queued, *, queue, target, common):
    c = common
    active = [s for s in state['submissions'] if s['status'] == 'submitted'
              and own_rows(queued, s['job_id'])]
    c.need(len(active) <= 1, 'More than one active extension array')
    if not active:
        return
    submission = active[0]
    job = submission['job_id']
    fresh = queue()
    own = own_rows(fresh, job)
    pending = state.get('throttle_intent')
    if pending:
        c.need(pending['array_job'] == job, 'Unresolved throttle belongs to another array')
        if pending.get('returncode') not in (None, 0):
            c.need(finished_members_only(pending.get('stderr', ''), job),
                   'Unknown previous throttle error requires review')
    if not any(r['state'] in WAITING for r in own):
        decision = dict(at=c.utc(), array_job=job, scheduling_only=True,
            action='skipped_no_pending_members', live_members=len(own),
            successful_throttle_update_claimed=False, scientific_parameters_unchanged=True)
        if pending:
            decision['superseded_unresolved_intent'] = pending
            observed = readback(submission)
            decision['readback'] = observed
            caps = {r['throttle'] for r in observed}
            c.need(len(caps) == 1, 'Inconsistent observed array throttles')
            submission['current_throttle'] = next(iter(caps))
            state.pop('throttle_intent')
        state['last_capacity_decision'] = decision
        state.setdefault('throttle_decisions', []).append(decision)
        c.write(control / 'state.json', state)
        return
    decision = target(fresh, job)
    current = submission.get('current_throttle', submission['concurrency'])
    decision.update(at=c.utc(), array_job=job, previous=current, scheduling_only=True,
        scientific_parameters_unchanged=True, unrelated_controllers_not_constrained=True)
    if decision['target'] != current or pending:
        if pending:
            decision['reconciles_unresolved_intent'] = pending
        argv = ['scontrol', 'update', 'JobId=' + job,
                'ArrayTaskThrottle=' + str(decision['target'])]
        decision['argv'] = argv
        state['throttle_intent'] = decision
        c.write(control / 'state.json', state)
        result = subprocess.run(argv, text=True, capture_output=True, timeout=45)
        decision.update(returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
        state.setdefault('throttle_decisions', []).append(decision)
        c.write(control / 'state.json', state)
        c.need(result.returncode == 0 or finished_members_only(result.stderr, job),
               'Own-array throttle update has an unrecognized failure')
        after = own_rows(queue(), job)
        observed = readback(submission)
        decision['readback'] = observed
        active_caps = {r['throttle'] for r in observed if r['state'] in WAITING | {'RUNNING','CONFIGURING','COMPLETING','SUSPENDED'}}
        if any(r['state'] in WAITING for r in after):
            c.need(active_caps == {decision['target']}, 'Live array throttle differs from requested target')
            submission['current_throttle'] = decision['target']
            decision['resolution'] = 'readback_verified_target'
        else:
            caps = {r['throttle'] for r in observed}
            c.need(len(caps) == 1, 'Inconsistent observed array throttles')
            submission['current_throttle'] = next(iter(caps))
            decision['resolution'] = 'superseded_no_pending_members'
        decision['successful_throttle_update_claimed'] = decision['resolution'] == 'readback_verified_target'
        state.pop('throttle_intent', None)
    state['last_capacity_decision'] = decision
    c.write(control / 'state.json', state)
