#!/usr/bin/env python3
"""Advance completed sample fits to evaluation/geometry checks; no scientific work here."""
import json
import subprocess
import time
from common import ROOT, OUT, samples, write_json, utc

ledger=OUT/'dispatch_ledger.json'
state=json.loads(ledger.read_text()) if ledger.exists() else {'evaluation':{},'quality':{}}
while True:
    queue=subprocess.run(['squeue','-r','-u','yimin','-h','-o','%i'],capture_output=True,text=True,check=True).stdout.splitlines()
    budget=max(0,980-len(queue))
    done={'evaluation':0,'quality':0}
    for sample in samples():
        prep=OUT/'prepared'/sample/'manifest.json'
        if not prep.exists():continue
        pm=json.loads(prep.read_text())
        if pm.get('status')!='completed':continue
        terminal=True
        for feature in ['all','hvg500','hvg1000','hvg2000','hvg3000','hvg5000','hvg2000_markers','seurat2000','seurat5000']:
            if feature in pm.get('unavailable_features',{}):continue
            path=OUT/'task_status'/sample/f'{feature}.json'
            if not path.exists() or json.loads(path.read_text()).get('status') not in ['completed','failed']:
                terminal=False;break
        for stage in ['evaluation','quality']:
            if (OUT/stage/sample/'COMPLETE').exists():
                done[stage]+=1;continue
            if not terminal or budget<1 or sample in state[stage]:continue
            script='evaluate' if stage=='evaluation' else 'quality'
            p=subprocess.run(['sbatch','--parsable',str(ROOT/f'handoff/hvg_ptc_20260916/{script}.sbatch'),
                              '--sample',sample],capture_output=True,text=True)
            if p.returncode:
                print(f'{utc()} submission delayed: {p.stderr.strip()}',flush=True)
                budget=0;break
            job=p.stdout.strip().split(';')[0]
            state[stage][sample]={'job':job,'submitted_at':utc()}
            budget-=1
            write_json(ledger,state)
            print(f'{utc()} {sample} {stage} {job}',flush=True)
    state.update(updated_at=utc(),completed=done)
    write_json(ledger,state)
    if done=={'evaluation':121,'quality':121}:
        print('All sample evaluations and geometry diagnostics complete.',flush=True)
        break
    time.sleep(30)
