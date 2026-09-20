"""Release bounded A1 SLURM waves only after documented pilot gates.

Administrative JSON/hash checks only; no matrices, scientific evaluation, or
outcome-dependent selection. Scientific runners remain frozen in embedding_v6.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
import socket
import subprocess
import time

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
CAMP = ROOT / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
OUT = CAMP / 'embedding'
HERE = OUT / 'dispatch'
SOURCE = CAMP / 'source_snapshots/embedding_v6'
SBATCH = CAMP / 'source_snapshots/embedding_workers4_profile_v1/job.sbatch'
PILOTS = ['TKU4163', 'NL022', 'SN040']
ACTIVE = {'PENDING', 'RUNNING', 'CONFIGURING', 'COMPLETING', 'REQUEUED', 'SUSPENDED'}


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f'.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')
    tmp.replace(path)


def checked(path, manifest='manifest.json', flag='COMPLETE'):
    try:
        return (path / flag).read_text().strip() == sha(path / manifest)
    except OSError:
        return False


def done(task):
    base = OUT / task['sample'] / task['budget'] / task['space']
    return (checked(base, 'fit_manifest.json', 'FIT_COMPLETE')
            and checked(base / 'evaluation') and checked(base / 'figures'))


def key(task):
    return '/'.join(task[k] for k in ('sample', 'budget', 'space'))


def release_gates():
    waiting, evidence = [], {}
    probe = CAMP / 'embedding_worker_probe/SN040/hvg2000/PCA2_workers4/validation.json'
    proof = json.loads(probe.read_text())
    assert proof['status'] == 'passed' and proof['operational_acceptance'] is True
    assert proof['DEG_and_density_R_objects_identical'] is True
    assert proof['numeric_scorer_sha256'] == sha(SOURCE / 'score_candidates.R')
    for path, digest in proof['source_manifests'].items():
        assert sha(path) == digest
    evidence[str(probe)] = sha(probe)
    # Six size/budget scopes exercise every frozen reducer and clusterer.
    for sample in PILOTS:
        for budget in ['hvg2000', 'hvg5000']:
            base = CAMP / 'embedding_geometry_smoke' / sample / budget
            if not checked(base):
                waiting.append(f'geometry smoke {sample}/{budget}')
                continue
            path = base / 'manifest.json'
            obj = json.loads(path.read_text())
            assert obj['status'] == 'passed' and obj['n_partitions'] == 91
            assert len(obj['spaces']) == 7 and all(s['n_conditions'] == 13 for s in obj['spaces'])
            assert obj['protocol_sha256'] == sha(CAMP / 'protocol/embedding.json')
            assert obj['numeric_source_sha256'] == sha(SOURCE / 'run.py')
            for record in obj['spaces']:
                unit = base / record['space']
                assert record['representation_sha256'] == sha(unit / 'representation.json')
                for route, digest in record['partition_manifests'].items():
                    name = route.removeprefix(record['space'] + '_')
                    partition = unit / name
                    assert digest == sha(partition / 'partition_manifest.json')
                    pm = json.loads((partition / 'partition_manifest.json').read_text())
                    assert pm['clusters_sha256'] == sha(partition / 'clusters.csv')
            evidence[str(path)] = sha(path)
    # Original-R scorer and terminal reproduction, independent of author truth.
    parity = OUT / 'TKU4163/hvg2000/anchor_parity/PARITY_PASSED.json'
    if not parity.exists():
        waiting.append('original anchor parity')
    else:
        obj = json.loads(parity.read_text())
        assert obj['status'] == 'passed'
        evidence[str(parity)] = sha(parity)
    # Independently audited end-to-end PCA2 at all three cell-count scales.
    for sample in PILOTS:
        task = dict(sample=sample, budget='hvg2000', space='PCA2')
        base = OUT / sample / 'hvg2000/PCA2'
        proof = CAMP / 'no_clustering/a1_first_review' / sample / 'hvg2000/PCA2/validation.json'
        if not done(task) or not proof.exists():
            waiting.append(f'end-to-end independent verification {sample}/PCA2')
            continue
        obj = json.loads(proof.read_text())
        assert obj['status'] == 'passed' and obj['n_conditions'] == 13 and obj['n_metric_rows'] == 39
        for field, path in [('fit_manifest_sha256', base / 'fit_manifest.json'),
                            ('evaluation_manifest_sha256', base / 'evaluation/manifest.json'),
                            ('figures_manifest_sha256', base / 'figures/manifest.json')]:
            assert obj[field] == sha(path), f'Stale pilot proof: {sample}/{field}'
        evidence[str(proof)] = sha(proof)
    return waiting, evidence


def run_command(args):
    return subprocess.run(args, capture_output=True, text=True, check=True, timeout=30).stdout


def scheduler(job_ids):
    if not job_ids:
        return {}
    output = run_command(['sacct', '--array', '-X', '-n', '-P', '-j', ','.join(job_ids),
                          '--format=JobID%80,State'])
    rows = {}
    for line in output.splitlines():
        jid, state, *_ = line.strip().split('|')
        rows[jid] = state.split()[0].rstrip('+')
    output = run_command(['squeue', '--array', '-h', '-j', ','.join(job_ids), '-o', '%i|%T'])
    for line in output.splitlines():
        jid, state = line.strip().split('|')
        rows[jid] = state
    return rows


def publish(state, message, waiting=(), error=None):
    total = len(state['tasks'])
    submitted = sum(len(w['tasks']) for w in state['waves'])
    value = dict(stage='A1_WAVES', work_package='A', status='needs_retry' if error else 'running',
        updated_at=utc(), summary=message,
        completed=[f'{submitted}/{total} remaining representation tasks submitted in bounded waves'],
        remaining=list(waiting), details=[
            'Each representation task evaluates all 13 frozen clustering candidates and produces every figure.',
            'At most two waves ×24 simultaneous jobs; new submissions keep all local-user queued/running jobs below900.',
            'Twenty-one original pilot representation units remain reserved to the existing serial pilot array.'],
        jobs=[dict(job_id=w['job_id'], purpose=f'A1 full wave {i}') for i, w in enumerate(state['waves']) if w.get('job_id')],
        evidence=[str(HERE / 'release.json')] if (HERE / 'release.json').exists() else [],
        whole_work_package_A_complete=False)
    if error:
        value['blockers'] = [error]
    write(HERE / 'status.json', value)


def initialize():
    path = HERE / 'state.json'
    if path.exists():
        return json.loads(path.read_text())
    spec = json.loads((CAMP / 'protocol/embedding.json').read_text())
    tasks = [dict(sample=s, budget=b, space=p) for s in spec['samples']
             for b in spec['budgets'] for p in spec['spaces']
             if not (s in PILOTS and b == 'hvg2000')]
    assert len(tasks) == 1673 and len({key(t) for t in tasks}) == 1673
    state = dict(created_at=utc(), tasks=tasks, waves=[],
        protocol_sha256=sha(CAMP / 'protocol/embedding.json'),
        frozen_source={str(p.relative_to(SOURCE)):sha(p) for p in SOURCE.rglob('*')
                       if p.is_file() and p.suffix in ('.py', '.R', '.sbatch')},
        sbatch_sha256=sha(SBATCH), max_wave_tasks=250, max_parallel_waves=2,
        concurrency_per_wave=24, max_total_submitted=900)
    write(path, state)
    return state


def iteration(state):
    assert state['protocol_sha256'] == sha(CAMP / 'protocol/embedding.json')
    assert state['sbatch_sha256'] == sha(SBATCH)
    for relative, digest in state['frozen_source'].items():
        assert sha(SOURCE / relative) == digest, 'Frozen numerical source was modified'
    assert not state.get('submission_intent'), 'Unresolved sbatch intent; inspect scheduler before any retry'
    waiting, evidence = release_gates()
    if waiting:
        publish(state, 'Full A1 waves waiting for independent pilot gates; existing jobs continue.', waiting)
        return False
    if not (HERE / 'release.json').exists():
        write(HERE / 'release.json', dict(status='passed', passed_at=utc(), evidence=evidence,
            scientific_protocol_sha256=state['protocol_sha256'],
            policy='All seven interfaces at three sizes/two budgets; all three PCA2 end-to-end independently verified.'))
    rows = scheduler([w['job_id'] for w in state['waves']])
    active = []
    for wave in state['waves']:
        taskfile = Path(wave['tasks_file'])
        assert sha(taskfile) == wave['tasks_sha256'], f'Changed wave task file: {taskfile}'
        assert json.loads(taskfile.read_text()) == wave['tasks'], f'Changed wave tasks: {taskfile}'
        own = [v for j,v in rows.items() if j == wave['job_id'] or j.startswith(wave['job_id'] + '_')]
        if all(done(t) for t in wave['tasks']):
            wave['status'] = 'artifacts_complete_pending_full_scientific_validation'
            if any(v in ACTIVE for v in own):
                active.append(wave)
            continue
        bad = [v for v in own if v not in ACTIVE | {'COMPLETED'}]
        assert not bad, f'Wave {wave["job_id"]} requires repair: {sorted(set(bad))}'
        if not own or all(v == 'COMPLETED' for v in own):
            raise RuntimeError(f'Wave {wave["job_id"]} has incomplete artifacts; refusing further submissions')
        active.append(wave)
    assigned = {key(t) for w in state['waves'] for t in w['tasks']}
    remaining = [t for t in state['tasks'] if key(t) not in assigned and not done(t)]
    if not remaining and not active:
        publish(state, 'All dispatched representation artifacts complete; full selection and independent final validation remain.')
        write(HERE / 'FIT_WAVES_COMPLETE.json', dict(status='completed', at=utc(), total=len(state['tasks'])))
        return True
    ready = [t for t in remaining if checked(OUT / t['sample'] / t['budget'] / 'geometry')]
    queue_count = len(run_command(['squeue', '--array', '-u', 'yimin', '-h', '-o', '%i']).splitlines())
    slots = max(0, min(state['max_wave_tasks'], state['max_total_submitted'] - queue_count))
    if len(active) < state['max_parallel_waves'] and ready and slots:
        wave_index = len(state['waves'])
        tasks = ready[:slots]
        taskpath = HERE / f'wave_{wave_index:02d}_tasks.json'
        write(taskpath, tasks)
        intent = dict(index=wave_index, tasks_file=str(taskpath), tasks_sha256=sha(taskpath), at=utc())
        state['submission_intent'] = intent
        write(HERE / 'state.json', state)
        job = run_command(['sbatch', '--parsable', f'--array=0-{len(tasks)-1}%24',
            '--time=1-00:00:00', f'--job-name=reviewA_wave{wave_index:02d}',
            f'--output={HERE}/wave_{wave_index:02d}_%A_%a.log', str(SBATCH), str(SOURCE), str(taskpath)]).strip().split(';')[0]
        assert job.isdigit(), f'Unexpected sbatch output: {job}'
        state['waves'].append(dict(job_id=job, tasks=tasks, **intent))
        del state['submission_intent']
    write(HERE / 'state.json', state)
    publish(state, 'Full A1 waves use frozen numerical code; completed outputs undergo final selection and verification.',
            [f'{len(remaining)} unassigned tasks at start of refresh; {len(ready)} have completed original-R geometry exports.'])
    return False


def main():
    assert os.environ.get('SLURM_JOB_ID'), 'Run administrative dispatcher within SLURM'
    parser = argparse.ArgumentParser();parser.add_argument('--once', action='store_true')
    args = parser.parse_args()
    HERE.mkdir(parents=True, exist_ok=True)
    # mkdir is atomic on GPFS across nodes; flock is not relied upon here.
    # A lock left by an abrupt job death requires an explicit owner audit.
    lock = HERE / 'manager.lock'
    lock.mkdir()
    write(lock / 'owner.json', dict(job=os.environ['SLURM_JOB_ID'],
        step=os.environ.get('SLURM_STEP_ID'), host=socket.gethostname(), pid=os.getpid(), started_at=utc()))
    try:
        state = initialize()
        while not (HERE / 'STOP').exists():
            try:
                finished = iteration(state)
            except Exception as exc:
                publish(state, 'A1 submission paused for a source/gate/job integrity issue; active scientific jobs are preserved.', error=repr(exc))
                raise
            if finished or args.once:
                break
            time.sleep(45)
    finally:
        (lock / 'owner.json').unlink()
        lock.rmdir()


if __name__ == '__main__':
    main()
