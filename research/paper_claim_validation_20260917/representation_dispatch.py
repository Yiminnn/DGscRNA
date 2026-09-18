"""Schedule only prespecified representation controls after default parity."""
import json
import time
from common import OUT, checked, write_json, utc
from dispatch import freeze, states, submit
from representation_controls import configurations

def tick():
    path=OUT/'representation_dispatch_state.json'
    state=json.loads(path.read_text()) if path.exists() else dict(jobs={})
    if not (OUT/'verification/representation_default_parity.json').exists():return False
    assert json.loads((OUT/'verification/representation_default_parity.json').read_text())['status']=='passed'
    active=states();source=freeze('representation_controls')
    slots=sum(v[0]=='claim_GBM_representation' for v in active.values())
    taskdir=OUT/'protocol/representation_controls';taskdir.mkdir(exist_ok=True)
    cases=[cfg for sample in (OUT/'protocol/pilot_samples.txt').read_text().split() for budget in ['hvg2000','hvg5000'] for cfg in configurations(sample,budget)]
    assert len(cases)==180
    for cfg in cases:
        key=cfg['sample']+'/'+cfg['budget']+'/'+cfg['name'];r=state['jobs'].setdefault(key,{})
        dest=OUT/'GBM_representation_controls'/key
        if checked(dest):r['status']='complete';continue
        if not checked(OUT/'GBM'/cfg['sample']/cfg['budget'],'prepare_manifest.json','PREPARED'):continue
        if r.get('job') in active:continue
        if r.get('job') and time.time()-r['submitted']<180:continue
        if r.get('job'):
            import subprocess
            raw=subprocess.check_output(['sacct','-j',r['job'],'-X','-n','-P','--format=JobID,State,ExitCode'],text=True)
            rows=[v.split('|') for v in raw.splitlines() if v.startswith(r['job']+'|')]
            if rows and rows[0][1].split()[0] in ['FAILED','OUT_OF_MEMORY','TIMEOUT','CANCELLED','COMPLETED','NODE_FAIL']:
                r.update(status='needs_review',accounting=rows[0])
            continue
        if slots>=24:continue
        task=taskdir/(key.replace('/','__')+'.json')
        if task.exists():assert json.loads(task.read_text())==cfg
        else:write_json(task,cfg)
        jid=submit(source,'representation_controls.py',[task],['--job-name=claim_GBM_representation','--cpus-per-task=4','--mem=40G','--time=04:00:00'])
        r.update(job=jid,submitted=time.time(),status='submitted',source=str(source));slots+=1;write_json(path,state)
    n=sum(r.get('status')=='complete' for r in state['jobs'].values())
    state.update(expected=180,completed=n,last_check=utc(),review_required={k:r['accounting'] for k,r in state['jobs'].items() if r.get('status')=='needs_review'})
    write_json(path,state);print('REPRESENTATION_PROGRESS',n,'/180',flush=True)
    return n==180

if __name__=='__main__':
    while not tick():time.sleep(45)
