"""Metadata-only contracts for the separately frozen GBM extension campaign."""
from pathlib import Path
import importlib.util
import os
import re

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
spec = importlib.util.spec_from_file_location('_frozen_core_metadata', PARENT / 'manage_core.py')
core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)
need, sha, load, write, utc = core.need, core.sha, core.load, core.write, core.utc
LIMITS = dict(queue_ceiling=900, queue_reserve=64, running_ceiling=256, core_running_reserve=16,
              initial_concurrency=16, minimum_concurrency=8, maximum_concurrency=64,
              max_own_outstanding=64, max_active_arrays=1, poll_seconds=60)
COUNTS = dict(geometry=363, mlp=54, neighbors=66, learning=18, no_cluster=242)
FULL_COUNTS = dict(COUNTS, representation=180, A1=1694)
PROFILES = dict(standard=dict(cpus=4, memory='64G', time='08:00:00'),
                A1=dict(cpus=4, memory='32G', time='24:00:00'))
OPTIONAL = dict(representation=dict(expected_tasks=180, enabled=False,
                    reason='Requires a frozen lock, full default parity, reviewed acceptance validator and new campaign gate'),
                A1=dict(expected_tasks=1694, enabled=False, resources=dict(cpus=4, memory='32G', time='24:00:00'),
                    reason='Requires full 91-terminal pilot, locked runtime, Lfine acceptance, and noDR-before-other-spaces dependencies'))
LAUNCHERS = ('common.py', 'receipts.py', 'worker.py', 'manage.py', 'array.sbatch')


def record(path):
    path = Path(path).resolve()
    return dict(path=str(path), sha256=sha(path))


def checked(record):
    path = Path(record['path']).resolve()
    need(path.is_file() and sha(path) == record['sha256'], f'Frozen file changed: {path}')
    return path


def receipt(directory, name, flag):
    path = Path(directory) / name
    need(path.is_file() and (Path(directory) / flag).read_text().strip() == sha(path),
         f'Incomplete or changed receipt: {path}')
    return load(path)


def lock_files(path):
    data = load(path)
    result = {}
    for key, item in data['files'].items():
        p, digest = (Path(item['path']), item['sha256']) if isinstance(item, dict) else (Path(key), item)
        need(p.is_absolute() and sha(p) == digest, f'Frozen adapter source changed: {p}')
        result[str(p.resolve())] = digest
    return result


def b_tasks():
    # Read the locked historical roster, without importing scientific modules.
    path = PARENT.parents[1] / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/controls/tasks.json'
    tasks = load(path)
    need(len(tasks) == 84 and sum(t['task'] == 'neighbors' for t in tasks) == 66
         and sum(t['task'] == 'learning' for t in tasks) == 18, 'Reviewer B roster differs')
    return tasks


def task_roster(release, core_tasks, *, full=False, representation_default=None):
    root = Path(release['output_root']).resolve()
    core_map = {(t['sample'], t['budget']): t for t in core_tasks}
    result = []

    def add(family, local, cfg, relative, script, extra=(), dependencies=None):
        sample, budget = cfg['sample'], cfg['budget']
        deps = dependencies or [(sample, budget)]
        output = root / 'extensions' / relative
        result.append(dict(index=len(result), family=family, local_index=local, sample=sample,
                           budget=budget, configuration=cfg, output=str(output),
                           script=str(PARENT / script), arguments=['--index', str(local), *extra],
                           dependencies=[core_map[key] for key in sorted(set(deps))],
                           extension_dependencies=[], profile='standard'))

    geometry = [t for t in core_tasks if t['budget'] in {'hvg2000', 'hvg5000', 'all'}]
    for i, t in enumerate(geometry):
        add('geometry', i, t, f"geometry_fixed_DL2000/{t['sample']}/{t['budget']}",
            'geometry_control.py', dependencies=[(t['sample'], t['budget']), (t['sample'], 'hvg2000')])
    mlp = load(PARENT / 'mlp_controls_tasks.json')
    for i, t in enumerate(mlp):
        need(t['index'] == i and t['n_terminal_conditions'] == 8, 'MLP task mismatch')
        relative = f"MLP_controls/{t['sample']}/{t['budget']}/{t['control']}"
        add('mlp', i, t, relative, 'mlp_controls_adapter.py',
            ['--out', str(root / 'extensions' / relative)])
    for family in ('neighbors', 'learning'):
        selected = [t for t in b_tasks() if t['task'] == family]
        for i, t in enumerate(selected):
            name = t['name'] if family == 'neighbors' else f"seed{t['model_seed']}"
            add(family, i, t, f"reviewer_b/{family}/{t['sample']}/{t['budget']}/{name}",
                f'reviewer_b_adapter/{family}.py')
    a2 = [t for t in core_tasks if t['budget'] in {'hvg2000', 'hvg5000'}]
    for i, t in enumerate(a2):
        add('no_cluster', i, t, f"no_cluster/{t['sample']}/{t['budget']}", 'no_cluster_adapter/run.py')
    if full:
        need(representation_default, 'Full campaign needs a verified representation default proof')
        representation = load(PARENT / 'representation_tasks.json')
        need(len(representation) == 180, 'Representation roster incomplete')
        for i, t in enumerate(representation):
            cfg = t['configuration']
            relative = f"representation_controls/{cfg['sample']}/{cfg['budget']}/{cfg['name']}"
            add('representation', i, cfg, relative, 'representation_adapter.py',
                ['--mode', 'task', '--default-proof', str(representation_default), '--out', str(root / 'extensions' / relative)])
        a1 = load(PARENT / 'embedding_adapter/tasks.json')
        need(len(a1) == 1694, 'A1 roster incomplete')
        offset = len(result)
        for i, t in enumerate(a1):
            need(t['index'] == i and t['n_partitions'] == 13, 'A1 task identity differs')
            add('A1', i, t, f"A1/{t['sample']}/{t['budget']}/{t['space']}", 'embedding_adapter/run.py')
            task = result[-1]
            task.update(profile='A1', arguments=[],
                        verification_output=str(root / 'extensions/A1_verification' / t['sample'] / t['budget'] / t['space']),
                        extension_dependencies=[] if t['space'] == 'noDR' else [offset + (i // 7) * 7])
            need(a1[(i // 7) * 7]['space'] == 'noDR' and
                 a1[(i // 7) * 7]['sample'] == t['sample'] and a1[(i // 7) * 7]['budget'] == t['budget'],
                 'A1 noDR export dependency differs')
    counts = FULL_COUNTS if full else COUNTS
    need({f: sum(t['family'] == f for t in result) for f in counts} == counts, 'Extension roster differs')
    need(len({t['output'] for t in result}) == len(result), 'Overlapping extension outputs')
    return result


def gate_metadata(path, *, candidate=False):
    gate = load(path)
    allowed = {'reviewed_extension_campaign'} | ({'candidate_pending_independent_review'} if candidate else set())
    need(gate.get('status') in allowed and gate.get('schema') == 1, 'Extension campaign gate is not reviewed')
    release_path = checked(gate['release_gate'])
    release, core_tasks, root = core.gate_metadata(release_path)
    full = gate['scope'] == 'full_2617'
    need(gate['scope'] in {'initial_743', 'full_2617'}, 'Unsupported campaign scope')
    need(gate['limits'] == LIMITS and gate['families'] == (FULL_COUNTS if full else COUNTS)
         and gate['optional_families'] == {name:dict(value, enabled=full) for name,value in OPTIONAL.items()}
         and gate['resource_profiles'] == PROFILES,
         'Campaign scope or scheduler limits changed')
    need(Path(gate['output_root']).resolve() == root, 'Campaign root differs from release')
    need(set(gate['launchers']) == set(LAUNCHERS), 'Missing extension launchers')
    for name, item in gate['launchers'].items():
        need(checked(item) == HERE / name, 'Wrong launcher path')
    for path_text, digest in gate['source_files'].items():
        checked(dict(path=path_text, sha256=digest))
    for item in gate['source_locks']:
        lock_files(checked(item))
    tasks = load(checked(gate['task_manifest']))
    default = gate.get('representation_default_proof')
    need(tasks == task_roster(release, core_tasks, full=full, representation_default=default), 'Task manifest differs from frozen scientific rosters')
    need(gate['task_count'] == len(tasks) == (2617 if full else 743), 'Incomplete extension campaign')
    for item in gate['pilots']:
        checked(item)
        proof = load(item['path'])
        wanted = 'passed' if item['role'] == 'representation_default' else 'passed_exact'
        need(proof.get('status') == wanted, 'A required full pilot did not pass its declared parity contract')
        for field, expected in item['expected'].items():
            need(proof.get(field) == expected, f'Pilot coverage differs: {item["path"]}: {field}')
    roles = {'geometry', 'mlp_original', 'mlp_seed0', 'B_neighbors', 'B_learning', 'A2'}
    if full:
        roles |= {'representation_default', 'representation_changed', 'A1_full91'}
        need(gate['A1_runtime'] == load(PARENT / 'embedding_adapter/RUNTIME.json'), 'A1 runtime lock differs')
        checked(gate['A1_runtime']['wheel'])
    need({p['role'] for p in gate['pilots']} == roles, 'Missing full pilots')
    return gate, tasks, release, root


def dependency_proofs(task, release, root, cache=None):
    cache = {} if cache is None else cache
    proofs = []
    for dep in task['dependencies']:
        key = dep['index']
        if key not in cache:
            cache[key] = core.acceptance(root, dep, sha(PARENT / 'RELEASE_GATE.json'), release['wheel']['sha256'])
        if cache[key] is None:
            return None
        proofs.append(cache[key])
    return proofs


def worker_receipt_root(root, index):
    return root / 'control/extension_manager/units' / f'{index:04d}'


def artifact_stats(paths):
    result = {}
    for path in paths:
        stat = Path(path).stat()
        result[path] = [stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns, stat.st_ino]
    return result


def historical_A1_ready(task, *, deep=False):
    """Return old-source receipt hashes only after every fitting endpoint exists."""
    if task['family'] != 'A1':
        return {}
    old = Path(task['configuration']['expected_reference'])
    if not (old / 'FIT_COMPLETE').exists():
        return None
    fit = receipt(old, 'fit_manifest.json', 'FIT_COMPLETE')
    cfg = fit['config']
    need(fit['status'] == 'fit_complete_evaluation_pending' and all(cfg[k] == task[k] for k in ['sample', 'budget'])
         and cfg['space'] == task['configuration']['space'] and len(cfg['conditions']) == 13
         and cfg['reference_labels_used_for_fit'] is False, 'Old A1 reference contract differs')
    receipt(old, 'representation.json', 'REPRESENTATION_COMPLETE')
    receipt(old.parent / 'geometry', 'manifest.json', 'COMPLETE')
    paths = [old / 'fit_manifest.json', old / 'representation.json', old.parent / 'geometry/manifest.json']
    for condition in cfg['conditions']:
        directory = Path(condition['dest']).resolve()
        need(directory.is_relative_to(old.resolve()), 'Old A1 condition escapes its space')
        partition = receipt(directory, 'partition_manifest.json', 'PARTITION_COMPLETE')
        score = receipt(directory, 'score_manifest.json', 'SCORE_COMPLETE')
        terminal = directory / 'terminal/L00_mean'
        tm = receipt(terminal, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
        need(tm['status'] == 'completed' and tm['terminal_valid'] is True, 'Old A1 terminal invalid')
        if deep:
            need(sha(directory / 'clusters.csv') == partition['clusters_sha256']
                 and sha(directory / 'initial_calls.csv.gz') == score['initial_sha256']
                 and sha(terminal / 'terminal.npz') == tm['terminal_sha256']
                 and sha(terminal / 'predictions.csv.gz') == tm['predictions_sha256'], 'Old A1 fitting payload changed')
        paths += [directory / 'partition_manifest.json', directory / 'score_manifest.json', terminal / 'terminal_manifest.json']
    return {str(p): sha(p) for p in paths}


def safe_identity():
    import pwd
    need(pwd.getpwuid(os.getuid()).pw_name == 'yimin', 'Only local yimin is authorized')
    need(os.environ.get('SLURM_JOB_ID'), 'Run within SLURM')
    need(not os.sys.flags.optimize, 'Python -O is prohibited')
