"""Parallel SCINA library tasks, each gated by its own three-size pilots."""
import csv
import json
import subprocess
import time
from common import OUT,checked,write_json,utc
from dispatch import freeze,submit,states

def tick():
    path=OUT/'SCINA_library_dispatch_state.json'
    state=json.loads(path.read_text()) if path.exists() else dict(jobs={})
    active=states();source=freeze('SCINA_libraries')
    slots=sum(name=='claim_GBM_SCINA_library' for name,status in active.values())
    core=json.loads((OUT/'dispatch_state.json').read_text())
    desired=96 if core.get('n_evaluated_plotted')==726 else 64 if core.get('n_terminal_complete',0)>=600 else 32
    # Leave a cushion for independent controllers; SLURM still enforces MaxJobsPU.
    other=sum(name!='claim_GBM_SCINA_library' for name,status in active.values())
    limit=max(0,min(desired,220-other))
    pilots=['TKU4163','NL022','SN040']
    gates={i:all(checked(OUT/'comparators/SCINA'/s/f'L{i:02d}') for s in pilots) for i in range(16)}
    rows=list(csv.DictReader((OUT/'protocol/sample_order.csv').open()))
    aux=json.loads((OUT/'aux_dispatch_state.json').read_text())
    for row in rows:
        sample=row['sample']
        # Preserve any whole-sample SCINA task submitted by an earlier scheduler.
        previous=aux['jobs'].get('SCINA/'+sample,{})
        if previous.get('job') in active:continue
        if not checked(OUT/'GBM'/sample/'hvg2000','prepare_manifest.json','PREPARED'):continue
        for i in range(16):
            key=f'{sample}/L{i:02d}';r=state['jobs'].setdefault(key,{})
            if checked(OUT/'comparators/SCINA'/key):r['status']='complete';continue
            if not gates[i] or r.get('job') in active:continue
            if r.get('job') and time.time()-r['submitted']<180:continue
            if r.get('job'):
                lines=subprocess.check_output(['sacct','-j',r['job'],'-X','-n','-P','--format=JobID,State,ExitCode'],text=True).splitlines()
                found=[v.split('|') for v in lines if v.startswith(r['job']+'|')]
                if found and found[0][1].split()[0] in ['FAILED','OUT_OF_MEMORY','TIMEOUT','CANCELLED','COMPLETED','NODE_FAIL']:
                    r.update(status='needs_review',accounting=found[0])
                continue
            if slots>=limit:continue
            job=submit(source,'scina_R.R',[sample,i],['--job-name=claim_GBM_SCINA_library','--cpus-per-task=2','--mem=32G','--time=04:00:00'])
            r.update(job=job,submitted=time.time(),status='submitted',source=str(source));slots+=1;write_json(path,state)
    completed=sum(r.get('status')=='complete' for r in state['jobs'].values())
    state.update(expected=121*16,completed=completed,concurrency_limit=limit,three_size_gate_by_library=gates,last_check=utc(),
                 review_required={k:r['accounting'] for k,r in state['jobs'].items() if r.get('status')=='needs_review'})
    write_json(path,state);print('SCINA_LIBRARY_PROGRESS',completed,'/1936',flush=True)
    return completed==121*16

if __name__=='__main__':
    while not tick():time.sleep(45)
