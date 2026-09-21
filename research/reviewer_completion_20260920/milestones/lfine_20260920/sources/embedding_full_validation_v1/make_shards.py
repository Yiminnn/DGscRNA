"""Build bounded immutable validation shards and print root submission commands."""
from pathlib import Path
import argparse,hashlib,json
import common as c
p=argparse.ArgumentParser();p.add_argument('--wave-tasks');p.add_argument('--pilot',action='store_true');p.add_argument('--after-job');p.add_argument('--concurrency',type=int,default=8)
a=p.parse_args();contract=c.contract();alltasks=c.js(c.OUT/'tasks.json')
assert 1<=a.concurrency<=16
assert not(a.wave_tasks and a.pilot)
if a.wave_tasks:
    selected=c.js(a.wave_tasks)
elif a.pilot:selected=c.pilot_tasks()
else:selected=alltasks
c.validate_tasks(selected)
if not a.pilot:c.pilot_accepted()
groups={}
for t in selected:
    if c.valid_report(t):continue
    groups.setdefault((t['sample'],t['budget']),[]).append(t)
shards=[dict(sample=k[0],budget=k[1],tasks=v) for k,v in groups.items()]
assert all(1<=len(s['tasks'])<=7 for s in shards)
identity=dict(shards=shards,contract_sha256=c.sha(c.OUT/'contract.json'),pilot=a.pilot,
    wave_tasks=str(Path(a.wave_tasks).resolve()) if a.wave_tasks else None,wave_tasks_sha256=c.sha(a.wave_tasks) if a.wave_tasks else None)
fingerprint=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:20]
bundle=c.OUT/'launches'/fingerprint;bundle.mkdir(parents=True,exist_ok=True)
path=bundle/'shards.json'
if path.exists():assert c.js(path)==identity
else:c.write(path,identity)
if a.after_job:assert a.after_job.isdigit()
cmd=['sbatch',f'--array=0-{len(shards)-1}%{a.concurrency}'] if shards else []
if cmd and a.after_job:cmd.append('--dependency=afterok:'+a.after_job)
if cmd:cmd += [f'--output={bundle}/%A_%a.log',str(c.CODE/'job.sbatch'),str(path)]
print(json.dumps(dict(shard_file=str(path),n_shards=len(shards),n_representations=sum(len(s['tasks']) for s in shards),
    command=cmd,submitted=False,note='Root enforces combined validation concurrency across simultaneously submitted bundles; no sbatch is executed here'),indent=2))
