"""Create an explicitly versioned scheduler fix without changing frozen workers."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent


def main():
    source = HERE / 'manage.py'
    text = source.read_text()
    replacements = []
    def replace(old, new):
        nonlocal text
        assert text.count(old) == 1
        text = text.replace(old, new)
        replacements.append(dict(old=old, new=new))
    start = text.index('def adapt_throttle(')
    end = text.index('\n\ndef recover_intent(', start)
    replace(text[start:end], '''def adapt_throttle(state, control, queued):
    import throttle_v2
    return throttle_v2.reconcile(state, control, queued, queue=queue,
                                target=throttle_target, common=c)


def controller_contract():
    path = Path(__file__).with_name('CONTROLLER_V2_LOCK.json')
    value = c.load(path)
    c.need(value['status'] == 'reviewed_scheduler_only_correction', 'Scheduler correction is not reviewed')
    for source, digest in value['files'].items():
        c.need(c.sha(source) == digest, 'Scheduler correction source changed: ' + source)
    c.need(value['previous_controller']['sha256'] == c.sha(Path(__file__).with_name('manage.py')),
           'Frozen original scheduler changed')
    return dict(path=str(path), sha256=c.sha(path))


def migrate_scheduler_state(state, control, gate_hash, tasks, batch_size):
    contract = controller_contract()
    if state['controller_sha256'] == c.sha(__file__):
        c.need(state['controller_override']['contract'] == contract, 'Scheduler correction contract changed')
        return
    import manage as original
    original.registry(state, gate_hash, tasks, batch_size)
    prior = control / 'controller_v2_migration_prior.json'
    if prior.exists():
        c.need(c.load(prior) == state, 'Inconsistent interrupted scheduler migration requires review')
    else:
        c.write(prior, state)
    state['controller_sha256'] = c.sha(__file__)
    state['controller_override'] = dict(contract=contract, previous_state=c.record(prior),
        migrated_at=c.utc(), reason='Read back applied throttle when completed array members produce scontrol nonzero',
        submissions_and_acceptances_preserved=True, campaign_gate_unchanged=True, scientific_workers_unchanged=True)
    c.write(control / 'state.json', state)
''')
    replace("    c.need(c.sha(gate_path) == gate_hash, 'Campaign gate changed during execution')",
        "    contract = controller_contract()\n    c.need(state['controller_override']['contract'] == contract, 'Scheduler correction contract changed')\n    c.need(c.sha(gate_path) == gate_hash, 'Campaign gate changed during execution')")
    replace("    outstanding = len(assigned - {int(k) for k in state['accepted']})",
        "    if not busy and state.get('throttle_intent'):\n        import throttle_v2\n        throttle_v2.retire_completed_intent(state, control, records, queued, c)\n    outstanding = len(assigned - {int(k) for k in state['accepted']})")
    replace("    args = parser.parse_args(); c.safe_identity()", "    args = parser.parse_args(); c.safe_identity()\n    contract = controller_contract()")
    replace("            state = c.load(control / 'state.json'); registry(state, gate_hash, tasks, args.batch_size)",
        "            state = c.load(control / 'state.json')\n            migrate_scheduler_state(state, control, gate_hash, tasks, args.batch_size)\n            registry(state, gate_hash, tasks, args.batch_size)")
    replace("                         limits=c.LIMITS, task_count=len(tasks), batch_size=args.batch_size, submissions=[], accepted={})",
        "                         limits=c.LIMITS, task_count=len(tasks), batch_size=args.batch_size, submissions=[], accepted={},\n                         controller_override=dict(contract=contract, fresh_controller=True))")
    replace("                         acceptances=state['accepted'], optional_families=gate['optional_families'], PTC_started=False)",
        "                         acceptances=state['accepted'], optional_families=gate['optional_families'], PTC_started=False,\n                         controller_override=state['controller_override'])")
    target = HERE / 'manage_v2.py'
    with target.open('x') as stream:
        stream.write(text)
    inverse = text
    for item in reversed(replacements):
        assert inverse.count(item['new']) == 1
        inverse = inverse.replace(item['new'], item['old'])
    assert inverse == source.read_text()
    manifest = dict(status='source_prepared_pending_regression_review', original=str(source),
        original_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        derivative=str(target), derivative_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        literal_replacements=replacements, exact_reverse_reconstruction=True,
        scientific_workers_or_campaign_gate_changed=False)
    with (HERE / 'CONTROLLER_V2_DERIVATION.json').open('x') as stream:
        json.dump(manifest, stream, indent=2); stream.write('\n')


if __name__ == '__main__':
    main()
