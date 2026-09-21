"""One sample/budget shard per SLURM task; at most seven units sequentially."""
from pathlib import Path
import hashlib,json,os,sys
import common as c
c.require_slurm();contract=c.contract();bundle=c.js(sys.argv[1])
assert bundle['contract_sha256']==c.sha(c.OUT/'contract.json')
assert Path(sys.argv[1]).resolve().parent.name==hashlib.sha256(json.dumps(bundle,sort_keys=True).encode()).hexdigest()[:20]
alltasks=[t for s in bundle['shards'] for t in s['tasks']];c.validate_tasks(alltasks)
assert all(len(s['tasks'])<=7 and all(t['sample']==s['sample'] and t['budget']==s['budget'] for t in s['tasks']) for s in bundle['shards'])
if bundle['pilot']:assert set(map(c.taskkey,alltasks))<=set(map(c.taskkey,c.pilot_tasks()))
else:c.pilot_accepted()
if bundle['wave_tasks']:
    assert c.sha(bundle['wave_tasks'])==bundle['wave_tasks_sha256']
    assert set(map(c.taskkey,alltasks))<=set(map(c.taskkey,c.js(bundle['wave_tasks'])))
shard=bundle['shards'][int(os.environ['SLURM_ARRAY_TASK_ID'])]
assert 1<=len(shard['tasks'])<=7
from validate_unit import run
for task in shard['tasks']:
    dest=c.reportpath(task);dest.mkdir(parents=True,exist_ok=True)
    if c.valid_report(task,deep=True):print('VERIFIED_UNIT_REUSED',task,flush=True);continue
    lock=dest/'ACTIVE';lock.mkdir()  # Fail closed on duplicate or stale owner.
    c.write(lock/'owner.json',dict(job=os.environ['SLURM_JOB_ID'],pid=os.getpid(),task=task))
    try:
        result=run(task,contract)
        print('INDEPENDENT_A1_UNIT_PASSED',task,flush=True)
    except Exception as e:
        c.write(dest/('FAILURE_'+os.environ['SLURM_JOB_ID']+'.json'),dict(status='failed',task=task,error=repr(e),at=c.utc(),model_fits=0))
        raise
    finally:
        (lock/'owner.json').unlink();lock.rmdir()
print('A1_VALIDATION_SHARD_COMPLETE',shard['sample'],shard['budget'],flush=True)
