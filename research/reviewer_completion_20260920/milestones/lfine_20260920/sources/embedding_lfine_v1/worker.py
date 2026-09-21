"""Evaluate and independently verify at most seven already accepted original units."""
from pathlib import Path
import hashlib,json,os,sys
import common as c
c.require_slurm();contract=c.contract();c.archive_gate(contract)
path=Path(sys.argv[1]);bundle=c.js(path)
assert bundle['contract_sha256']==c.sha(c.OUT/'contract.json')
assert path.parent.name==hashlib.sha256(json.dumps(bundle,sort_keys=True).encode()).hexdigest()[:20]
shards=c.js(c.OUT/'shards.json');indices={(t['sample'],t['budget']):i for i,t in enumerate(shards)}
key=lambda t:(t['sample'],t['budget'],t['space'])
alltasks=[t for tasks in bundle['shards'] for t in tasks]
assert len({key(t) for t in alltasks})==len(alltasks)
for tasks in bundle['shards']:
    assert 1<=len(tasks)<=7 and len({(t['sample'],t['budget']) for t in tasks})==1
    for t in tasks:assert t['space'] in shards[indices[t['sample'],t['budget']]]['spaces']
if bundle['pilot']:assert {key(t) for t in alltasks}<={key(t) for t in contract['pilot_units']}
else:
    gate=c.js(c.OUT/'PILOT_ACCEPTANCE.json');assert gate['status']=='passed' and gate['contract_sha256']==c.sha(c.OUT/'contract.json')
    assert set(gate['proofs'])=={str(c.OUT/'verification'/t['sample']/t['budget']/t['space']/'manifest.json') for t in contract['pilot_units']}
    for path,digest in gate['proofs'].items():assert c.sha(path)==digest and c.checked_outputs(Path(path).parent)['status']=='passed'
from evaluate import run as evaluate
from verify_shard import run as verify
from selftest import run as selftest
assert selftest()['status']=='passed'
for task in bundle['shards'][int(os.environ['SLURM_ARRAY_TASK_ID'])]:
    c.unit_acceptance(task['sample'],task['budget'],task['space'],None)
    parent=c.OUT/'locks'/task['sample']/task['budget'];parent.mkdir(parents=True,exist_ok=True);lock=parent/(task['space']+'.ACTIVE');lock.mkdir()
    c.write(lock/'owner.json',dict(job=os.environ['SLURM_JOB_ID'],pid=os.getpid(),task=task))
    try:
        index=indices[task['sample'],task['budget']];evaluate(index,task['space']);verify(index,task['space'])
        print('A1_LFINE_UNIT_VERIFIED',task,flush=True)
    except Exception as error:
        c.write(parent/(task['space']+'.FAILURE.'+os.environ['SLURM_JOB_ID']+'.json'),dict(status='failed',error=repr(error),task=task,at=c.utc()));raise
    finally:(lock/'owner.json').unlink();lock.rmdir()
