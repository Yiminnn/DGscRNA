"""Metadata-only launch bundles; print commands, never invoke sbatch."""
from pathlib import Path
import argparse,hashlib,json
import common as c
p=argparse.ArgumentParser();p.add_argument('--pilot',action='store_true');p.add_argument('--tasks');p.add_argument('--concurrency',type=int,default=8);p.add_argument('--after-job');a=p.parse_args()
contract=c.contract();assert 1<=a.concurrency<=16 and not(a.pilot and a.tasks)
roster={(s['sample'],s['budget'],p) for s in c.js(c.OUT/'shards.json') for p in s['spaces']}
key=lambda t:(t['sample'],t['budget'],t['space'])
if a.pilot:tasks=contract['pilot_units']
elif a.tasks:tasks=c.js(a.tasks)
else:tasks=[dict(sample=s,budget=b,space=p) for s,b,p in sorted(roster)]
assert tasks and len({key(t) for t in tasks})==len(tasks) and all(set(t)=={'sample','budget','space'} and key(t) in roster for t in tasks)
if not a.pilot:
    gate=c.js(c.OUT/'PILOT_ACCEPTANCE.json');assert gate['status']=='passed' and gate['contract_sha256']==c.sha(c.OUT/'contract.json')
    expected={str(c.OUT/'verification'/t['sample']/t['budget']/t['space']/'manifest.json') for t in contract['pilot_units']}
    assert set(gate['proofs'])==expected
    for path,digest in gate['proofs'].items():assert c.sha(path)==digest and c.checked_outputs(Path(path).parent)['status']=='passed'
groups={}
for task in tasks:
    # Metadata proof checks are a prerequisite even for a delayed submission.
    c.unit_acceptance(*key(task),None)
    done=c.OUT/'verification'/task['sample']/task['budget']/task['space']
    if c.checked(done):c.checked_outputs(done);continue
    groups.setdefault((task['sample'],task['budget']),[]).append(task)
identity=dict(contract_sha256=c.sha(c.OUT/'contract.json'),pilot=a.pilot,shards=list(groups.values()))
name=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:20]
directory=c.OUT/'launches'/name;directory.mkdir(parents=True,exist_ok=True);path=directory/'tasks.json'
if path.exists():assert c.js(path)==identity
else:c.write(path,identity)
command=[]
if groups:
    command=['sbatch',f'--array=0-{len(groups)-1}%{a.concurrency}',f'--output={directory}/%A_%a.log']
    if a.after_job:assert a.after_job.isdigit();command.append('--dependency=afterok:'+a.after_job)
    command += [str(c.CODE/'job.sbatch'),str(path)]
print(json.dumps(dict(bundle=str(path),n_shards=len(groups),command=command,submitted=False),indent=2))
