"""Release bounded A1 SLURM waves only after documented pilot gates.

Administrative JSON/hash checks only; no matrices, scientific evaluation, or
outcome-dependent selection. Adaptive ICA is an explicit repaired comparator;
original parallel-solver failures and prior successful producers are retained.
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
SOURCE = CAMP / 'source_snapshots/embedding_v8'
SBATCH = CAMP / 'source_snapshots/embedding_workers4_profile_v2/job.sbatch'
PRESENTATION = CAMP / 'source_snapshots/embedding_presentation_v1'
REPAIR_POLICY = CAMP / 'protocol/embedding_convergence_repair_20260920_v2.json'
ACTIVATION = HERE / 'v3_activation.json'
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
    primary = (checked(base, 'fit_manifest.json', 'FIT_COMPLETE')
               and checked(base / 'evaluation') and checked(base / 'figures'))
    if not primary or task['space'] != 'ICA2':
        return primary
    if not checked(base / 'figures_adaptive_v1'):
        return False
    manifest = json.loads((base / 'figures_adaptive_v1/manifest.json').read_text())
    for field, path in [
        ('representation_manifest_sha256', base / 'representation.json'),
        ('fit_manifest_sha256', base / 'fit_manifest.json'),
        ('evaluation_manifest_sha256', base / 'evaluation/manifest.json'),
        ('policy_sha256', REPAIR_POLICY),
        ('presentation_source_manifest_sha256', PRESENTATION / 'SOURCE_MANIFEST.json'),
    ]:
        assert manifest[field] == sha(path), f'Stale adaptive presentation: {base}/{field}'
    assert manifest['space_display'] == 'ICA2 adaptive'
    assert manifest['every_candidate_plotted'] is True
    assert len(manifest['conditions']) == 13 and len(manifest['files']) == 26
    # Final selector independently verifies all 26 figure bytes per unit. This
    # scheduling completion check is not promoted to scientific acceptance.
    return True


def key(task):
    return '/'.join(task[k] for k in ('sample', 'budget', 'space'))


def release_gates():
    waiting, evidence = [], {}
    # Activation is written only after the independently reviewed frozen source,
    # successful-ICA regression, and completed-cache tamper tests all pass.
    if not ACTIVATION.exists():
        waiting.append('independent v8 source/cache/ICA acceptance')
    else:
        approval = json.loads(ACTIVATION.read_text())
        assert approval['status'] == 'passed'
        assert approval['source_manifest_sha256'] == sha(SOURCE / 'SOURCE_MANIFEST.json')
        assert approval['numeric_source_sha256'] == sha(SOURCE / 'run.py')
        assert approval['repair_policy_sha256'] == sha(REPAIR_POLICY)
        for name, digest in approval['evidence'].items():
            assert sha(name) == digest
        evidence[str(ACTIVATION)] = sha(ACTIVATION)
    presentation_proof = CAMP / 'embedding_v8_regression/presentation_validation.json'
    if not presentation_proof.exists():
        waiting.append('adaptive presentation validation')
    else:
        proof = json.loads(presentation_proof.read_text())
        assert proof['status'] == 'passed'
        assert proof['source_manifest_sha256'] == sha(PRESENTATION / 'SOURCE_MANIFEST.json')
        assert proof['n_figures'] == 26 and proof['all_original_figure_hashes_unchanged'] is True
        assert proof['figure_manifest_sha256'] == sha(Path(proof['unit']) / 'figures_adaptive_v1/manifest.json')
        evidence[str(presentation_proof)] = sha(presentation_proof)
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
            base = CAMP / 'embedding_geometry_smoke_v2' / sample / budget
            if not checked(base):
                waiting.append(f'geometry smoke {sample}/{budget}')
                continue
            path = base / 'manifest.json'
            obj = json.loads(path.read_text())
            assert obj['status'] == 'passed' and obj['n_partitions'] == 91
            assert len(obj['spaces']) == 7 and all(s['n_conditions'] == 13 for s in obj['spaces'])
            assert obj['protocol_sha256'] == sha(CAMP / 'protocol/embedding.json')
            assert obj['numeric_source_sha256'] == sha(SOURCE / 'run.py')
            assert obj['repair_policy_sha256'] == sha(REPAIR_POLICY)
            for record in obj['spaces']:
                unit = Path(record['unit_path'])
                assert unit.resolve().is_relative_to(base.resolve())
                assert record['representation_sha256'] == sha(unit / 'representation.json')
                representation = json.loads((unit / 'representation.json').read_text())
                if record['space'] == 'noDR':
                    assert representation['embedding_sha256'] is None
                else:
                    assert representation['embedding_sha256'] == sha(unit / 'embedding.csv')
                for route, digest in record['partition_manifests'].items():
                    name = route.removeprefix(record['space'] + '_')
                    partition = unit / name
                    assert digest == sha(partition / 'partition_manifest.json')
                    pm = json.loads((partition / 'partition_manifest.json').read_text())
                    assert pm['clusters_sha256'] == sha(partition / 'clusters.csv')
            evidence[str(path)] = sha(path)
    independent_smoke = CAMP / 'embedding_geometry_smoke_v2/gate_verification.json'
    if not independent_smoke.exists():
        waiting.append('independent aggregate smoke verification')
    else:
        proof = json.loads(independent_smoke.read_text())
        assert proof['status'] == 'passed' and proof['n_units'] == 6
        assert proof['n_spaces'] == 42 and proof['n_partitions'] == 546
        assert proof['numeric_source_sha256'] == sha(SOURCE / 'run.py')
        assert proof['source_manifest_sha256'] == sha(SOURCE / 'SOURCE_MANIFEST.json')
        assert proof['repair_policy_sha256'] == sha(REPAIR_POLICY)
        assert proof['every_cell_order_embedding_finite_partition_integer_verified'] is True
        assert proof['historical_files_unchanged'] is True
        for unit in proof['units']:
            assert unit['manifest_sha256'] == sha(Path(unit['unit_path']) / 'manifest.json')
        evidence[str(independent_smoke)] = sha(independent_smoke)
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


def publish(state, message, waiting=(), error=None, status='running'):
    total = len(state['tasks'])
    submitted = sum(len(w['tasks']) for w in state['waves'])
    value = dict(stage='A1_WAVES', work_package='A', status='needs_retry' if error else status,
        updated_at=utc(), summary=message,
        completed=[f'{submitted}/{total} remaining representation tasks submitted in bounded waves'],
        remaining=list(waiting), details=[
            'Each representation task evaluates all 13 frozen clustering candidates and produces every figure.',
            f"At most {state['max_parallel_waves']} waves ×{state['concurrency_per_wave']} simultaneous jobs; all local-user queued/running jobs total at most {state['max_total_submitted']} at submission.",
            'Twenty-one original pilot representation units remain reserved to their existing pilot and checkpoint-replacement jobs.',
            'ICA2 adaptive: parallel 5000/50000/100000, then deflation fallback; actual solver, cap and unsuccessful original attempts retained.'],
        jobs=[dict(job_id=w['job_id'], purpose=f'A1 full wave {i}') for i, w in enumerate(state['waves']) if w.get('job_id')],
        evidence=[str(HERE / 'release.json')] if (HERE / 'release.json').exists() else [],
        whole_work_package_A_complete=False)
    if error:
        value['blockers'] = [error]
    write(HERE / 'status.json', value)


def expected_tasks():
    spec = json.loads((CAMP / 'protocol/embedding.json').read_text())
    tasks = [dict(sample=s, budget=b, space=p) for s in spec['samples']
             for b in spec['budgets'] for p in spec['spaces']
             if not (s in PILOTS and b == 'hvg2000')]
    assert len(tasks) == 1673 and len({key(t) for t in tasks}) == 1673
    return tasks


def validate_state(state):
    assert state['tasks'] == expected_tasks(), 'Changed full task roster or ordering'
    assert (state['max_wave_tasks'], state['max_parallel_waves'],
            state['concurrency_per_wave'], state['max_total_submitted']) == (128, 6, 32, 900)
    permitted = {key(task) for task in state['tasks']}
    assigned = [key(task) for wave in state['waves'] for task in wave['tasks']]
    assert len(assigned) == len(set(assigned)), 'Overlapping wave task assignments'
    assert set(assigned) <= permitted, 'Unknown or reserved pilot task assigned'
    assert all(0 < len(wave['tasks']) <= 128 for wave in state['waves'])
    job_ids = [wave['job_id'] for wave in state['waves']]
    assert len(job_ids) == len(set(job_ids)), 'Repeated wave job identifier'


def initialize():
    path = HERE / 'state.json'
    if path.exists():
        state = json.loads(path.read_text())
        validate_state(state)
        return state
    tasks = expected_tasks()
    state = dict(created_at=utc(), tasks=tasks, waves=[],
        protocol_sha256=sha(CAMP / 'protocol/embedding.json'),
        frozen_source={str(p.relative_to(SOURCE)):sha(p) for p in SOURCE.rglob('*')
                       if p.is_file() and p.suffix in ('.py', '.R', '.sbatch')},
        source_manifest_sha256=sha(SOURCE / 'SOURCE_MANIFEST.json'),
        presentation_source_manifest_sha256=sha(PRESENTATION / 'SOURCE_MANIFEST.json'),
        repair_policy_sha256=sha(REPAIR_POLICY),
        sbatch_sha256=sha(SBATCH), max_wave_tasks=128, max_parallel_waves=6,
        concurrency_per_wave=32, max_total_submitted=900)
    write(path, state)
    validate_state(state)
    return state


def iteration(state):
    validate_state(state)
    assert state['protocol_sha256'] == sha(CAMP / 'protocol/embedding.json')
    assert state['sbatch_sha256'] == sha(SBATCH)
    assert state['source_manifest_sha256'] == sha(SOURCE / 'SOURCE_MANIFEST.json')
    assert state['presentation_source_manifest_sha256'] == sha(PRESENTATION / 'SOURCE_MANIFEST.json')
    assert state['repair_policy_sha256'] == sha(REPAIR_POLICY)
    for relative, digest in state['frozen_source'].items():
        assert sha(SOURCE / relative) == digest, 'Frozen numerical source was modified'
    for relative, digest in json.loads((SOURCE / 'SOURCE_MANIFEST.json').read_text()).items():
        assert sha(SOURCE / relative) == digest, 'Frozen numerical dependency was modified'
    for relative, digest in json.loads((PRESENTATION / 'SOURCE_MANIFEST.json').read_text()).items():
        assert sha(PRESENTATION / relative) == digest, 'Frozen presentation dependency was modified'
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
        publish(state, 'All dispatched representation artifacts complete; full selection and independent final validation remain.', status='awaiting_validation')
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
        job = run_command(['sbatch', '--parsable', f"--array=0-{len(tasks)-1}%{state['concurrency_per_wave']}",
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
