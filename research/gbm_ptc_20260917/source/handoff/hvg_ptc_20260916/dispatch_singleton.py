#!/usr/bin/env python3
"""Bounded dispatch of the verified secondary scoring sensitivity after primary evaluation."""
import json
import subprocess
import time
from common import ROOT,OUT,samples,write_json,utc

gate=OUT/'singleton_robustness_summary/manifest.json'
assert gate.exists() and json.loads(gate.read_text())['n_repairs_verified']>=115, 'Verified pilot required'
path=OUT/'singleton_dispatch_ledger.json'
state=json.loads(path.read_text()) if path.exists() else {'jobs':{}}
while True:
    lines=subprocess.run(['squeue','-r','-u','yimin','-h','-o','%i|%j'],capture_output=True,text=True,check=True).stdout.splitlines()
    active={line.split('|')[0] for line in lines}
    pending=sum(r['job'] in active for r in state['jobs'].values())
    budget=min(max(0,24-pending),max(0,980-len(lines)))
    complete=0
    for sample in samples():
        if (OUT/'singleton_robustness'/sample/'COMPLETE').exists():complete+=1;continue
        if budget<1 or sample in state['jobs'] or not (OUT/'evaluation'/sample/'COMPLETE').exists():continue
        p=subprocess.run(['sbatch','--parsable',str(ROOT/'handoff/hvg_ptc_20260916/singleton_robustness.sbatch'),
                          '--sample',sample],capture_output=True,text=True)
        if p.returncode:
            print(p.stderr,flush=True);budget=0;break
        job=p.stdout.strip().split(';')[0]
        state['jobs'][sample]={'job':job,'submitted_at':utc()};budget-=1
        write_json(path,state)
        print(utc(),sample,job,flush=True)
    state.update(updated_at=utc(),complete_samples=complete)
    write_json(path,state)
    if complete==121:break
    time.sleep(30)
