"""Build a nonlaunchable extension gate after validating locked pilot receipts."""
from pathlib import Path
import argparse
import common as c
import receipts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release-gate', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--full', action='store_true', help='Require all full pilots and build the unified2617 roster')
    parser.add_argument('--representation-default', type=Path)
    parser.add_argument('--representation-changed', type=Path)
    parser.add_argument('--A1-pilot', dest='a1_pilot', type=Path)
    args = parser.parse_args(); c.safe_identity()
    out = args.out.resolve()
    c.need(not out.exists(), 'Preserve existing candidate gate')
    release_path = args.release_gate.resolve()
    release, core_tasks, root = c.core.gate_metadata(release_path)
    if args.full:
        c.need(args.representation_default and args.representation_changed and args.a1_pilot, 'All optional full proofs required')
        args.representation_default = args.representation_default.resolve()
        args.representation_changed = args.representation_changed.resolve()
        args.a1_pilot = args.a1_pilot.resolve()
    tasks = c.task_roster(release, core_tasks, full=args.full, representation_default=args.representation_default)
    sources, locks = {}, []
    for relative in ['mlp_controls_sources.json', 'reviewer_b_adapter/SOURCE_LOCK.json', 'no_cluster_adapter/SOURCE_LOCK.json']:
        path = c.PARENT / relative
        locks.append(c.record(path)); sources.update(c.lock_files(path))
    sources[str(c.PARENT / 'geometry_control.py')] = c.sha(c.PARENT / 'geometry_control.py')
    sources[str(c.PARENT / 'mlp_controls_tasks.json')] = c.sha(c.PARENT / 'mlp_controls_tasks.json')
    b_roster = c.PARENT.parents[1] / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/controls/tasks.json'
    sources[str(b_roster)] = c.sha(b_roster)
    pilots = []

    def add_pilot(role, path, expected):
        item = dict(c.record(path), role=role, expected=expected)
        value = c.load(path)
        wanted = 'passed' if role == 'representation_default' else 'passed_exact'
        c.need(value['status'] == wanted and all(value.get(k) == v for k,v in expected.items()),
               f'Pilot failed/incomplete: {role}')
        pilots.append(item)

    checks = []
    for family, relative, selector in [
        ('geometry', 'pilots/geometry_fixed_DL2000_TKU4163_hvg2000', lambda t: t['family'] == 'geometry'),
        ('mlp', 'extensions/MLP_pilots/TKU4163_hvg2000_original', lambda t: t['family'] == 'mlp' and t['configuration']['control'] == 'original'),
        ('mlp', 'extensions/MLP_pilots/TKU4163_hvg2000_model_seed0', lambda t: t['family'] == 'mlp' and t['configuration']['control'] == 'model_seed0'),
        ('no_cluster', 'pilots/no_cluster_TKU4163_hvg2000', lambda t: t['family'] == 'no_cluster')]:
        task = next(t for t in tasks if t['sample'] == 'TKU4163' and t['budget'] == 'hvg2000' and selector(t))
        task = dict(task, output=str(root / relative))
        proof = receipts.validate(task, c.sha(release_path), deep=True, pilot=True)
        sources.update(proof['artifacts']); checks.append(proof)
        if family == 'geometry':
            add_pilot('geometry', root / relative / 'geometry_acceptance.json', dict(terminal_conditions=8, lfine_valid_threshold_rows=16))
        elif family == 'mlp':
            control = task['configuration']['control']
            add_pilot('mlp_original' if control == 'original' else 'mlp_seed0', root / relative / 'manifest.json',
                      dict(terminal_conditions=8, lfine_valid_threshold_rows=16, control=control))
        else:
            add_pilot('A2', root / relative / 'acceptance.json', dict(terminal_conditions=5, lfine_rows=10))
    bpilot = root / 'pilots/reviewer_b_default_TKU4163_hvg2000'
    for task in tasks:
        if task['family'] != 'neighbors' or task['sample'] != 'TKU4163' or task['budget'] != 'hvg2000':
            continue
        cfg = task['configuration']
        if cfg['snn_k'] != 20 or cfg['umap_neighbors'] != 30:
            continue
        pilot_task = dict(task, output=str(bpilot / 'neighbors/TKU4163/hvg2000' / cfg['name']))
        proof = receipts.validate(pilot_task, c.sha(release_path), deep=True, pilot=True)
        sources.update(proof['artifacts']); checks.append(proof)
    task = next(t for t in tasks if t['family'] == 'learning' and t['sample'] == 'TKU4163'
                and t['budget'] == 'hvg2000' and t['configuration']['model_seed'] == 42)
    task = dict(task, output=str(bpilot / 'learning/TKU4163/hvg2000/seed42'), pilot_epochs=10)
    proof = receipts.validate(task, c.sha(release_path), deep=True, pilot=True)
    sources.update(proof['artifacts']); checks.append(proof)
    add_pilot('B_neighbors', bpilot / 'neighbors_parity.json', dict(default_routes=3, terminal_conditions=6, Lfine_threshold_rows=12))
    add_pilot('B_learning', bpilot / 'learning_parity.json', dict(epochs=10))
    if args.full:
        for relative in ['representation_sources.json', 'embedding_adapter/ADAPTER_LOCK.json']:
            path = c.PARENT / relative
            locks.append(c.record(path)); sources.update(c.lock_files(path))
        sources[str(c.PARENT / 'representation_tasks.json')] = c.sha(c.PARENT / 'representation_tasks.json')
        for relative in ['embedding_adapter/tasks.json', 'embedding_adapter/RUNTIME.json']:
            sources[str(c.PARENT / relative)] = c.sha(c.PARENT / relative)
        default = c.receipt(args.representation_default.parent, 'manifest.json', 'COMPLETE')
        task = dict(index=-1, family='representation', sample=default['sample'], budget=default['budget'],
                    configuration={}, output=str(args.representation_default.parent), pilot_mode='default')
        proof = receipts.validate(task, c.sha(release_path), deep=True, pilot=True)
        sources.update(proof['artifacts']); checks.append(proof)
        add_pilot('representation_default', args.representation_default, dict(mode='default', terminal_conditions=8, lfine_valid_threshold_rows=16))
        changed = c.load(args.representation_changed)
        cfg = changed['configurations'][0]
        task = next(t for t in tasks if t['family'] == 'representation' and t['configuration'] == cfg)
        task = dict(task, output=str(args.representation_changed.parent))
        proof = receipts.validate(task, c.sha(release_path), deep=True, pilot=True)
        sources.update(proof['artifacts']); checks.append(proof)
        add_pilot('representation_changed', args.representation_changed, dict(mode='changed', terminal_conditions=2, lfine_valid_threshold_rows=4))
        a1 = c.receipt(args.a1_pilot.parent, 'manifest.json', 'COMPLETE')
        add_pilot('A1_full91', args.a1_pilot,
                  dict(n_representations=7, n_partitions=91, n_terminal_conditions=91, n_Lfine_threshold_rows=182))
        for task in [t for t in tasks if t['family'] == 'A1' and t['sample'] == a1['sample'] and t['budget'] == a1['budget']]:
            space = task['configuration']['space']
            specification = c.load(args.a1_pilot.parent / space / 'evaluation_spec.json')
            td = Path(specification['conditions'][0]['terminal_directory'])
            # The terminal directory is <space>/<partition>/terminal/<arm>.
            actual_space = td.parents[2]
            task = dict(task, output=str(actual_space), verification_output=str(args.a1_pilot.parent))
            proof = receipts.validate(task, c.sha(release_path), deep=True, pilot=True)
            sources.update(proof['artifacts']); checks.append(proof)
    # Pilot scientific payloads were checked above once. Keep only metadata
    # bindings in the frequently validated gate, avoiding repeated model IO.
    sources = {p:d for p,d in sources.items() if Path(p).suffix not in {'.npz', '.pt', '.gz', '.npy', '.bin'}}
    roster_path = out.with_name(out.stem + '.tasks.json')
    c.need(not roster_path.exists(), 'Preserve existing task manifest')
    c.write(roster_path, tasks)
    value = dict(schema=1, status='candidate_pending_independent_review', created_at=c.utc(),
                 release_gate=c.record(release_path), output_root=str(root), limits=c.LIMITS,
                 scope='full_2617' if args.full else 'initial_743',
                 families=c.FULL_COUNTS if args.full else c.COUNTS,
                 optional_families={name:dict(value, enabled=args.full) for name,value in c.OPTIONAL.items()}, task_count=len(tasks),
                 task_manifest=c.record(roster_path), launchers={name:c.record(c.HERE / name) for name in c.LAUNCHERS},
                 source_locks=locks, source_files=sources, pilots=pilots, pilot_validation=checks,
                 resource_profiles=c.PROFILES,
                 representation_default_proof=str(args.representation_default) if args.full else None,
                 A1_runtime=c.load(c.PARENT / 'embedding_adapter/RUNTIME.json') if args.full else None,
                 scientific_retry=False, automatic_optional_activation=False, automatic_PTC=False)
    c.write(out, value)
    c.gate_metadata(out, candidate=True)
    print(f'CANDIDATE_ONLY {len(tasks)} {out} {c.sha(out)}')


if __name__ == '__main__':
    main()
