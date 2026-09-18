"""Release larger cold-input resource runs only after measured smaller pilots."""
import json
import re
import shutil
import subprocess
import time
from common import OUT, checked, write_json, utc
from dispatch import freeze,submit,states

SIZES=[10000,30000,50000,100000,120000]
PILOTS={'DG-scRNA':'7363515','SCINA':'7363795','scDeepSort':'7363796'}

def dest(method,n):
    p=OUT/'scalability'/str(n)
    return p if method=='DG-scRNA' else p/method

def accounting(job):
    lines=subprocess.check_output(['sacct','-j',str(job),'-P','--format=JobID,State,ElapsedRaw,MaxRSS,ExitCode'],text=True).splitlines()
    data=[dict(zip(lines[0].split('|'),r.split('|'))) for r in lines[1:]]
    return data

def tick():
    path=OUT/'scalability_dispatch_state.json'
    state=json.loads(path.read_text()) if path.exists() else dict(jobs={})
    active=states();source=freeze('scalability')
    slots=sum((name.startswith('claim_scale_') or name.startswith('claim_scalability_')) and name!='claim_scale_dispatch' for name,status in active.values())
    for method in PILOTS:
        for n in SIZES:
            key=f'{method}/{n}';record=state['jobs'].setdefault(key,{})
            if n==10000 and not record:record.update(job=PILOTS[method],submitted=0,status='submitted',pilot=True)
            if checked(dest(method,n)):
                record['status']='complete'
                if 'accounting' not in record:record['accounting']=accounting(record['job'])
                continue
            if record.get('job') in active:continue
            if record.get('job') and time.time()-record.get('submitted',0)<180:continue
            if record.get('job'):
                acc=accounting(record['job'])
                exact=[r for r in acc if r['JobID']==record['job']]
                if exact and exact[0]['State'].split()[0] in ['FAILED','TIMEOUT','OUT_OF_MEMORY','CANCELLED','COMPLETED','NODE_FAIL']:
                    record.update(status='needs_review',accounting=acc)
                continue
            gate=10000 if n<=50000 else 50000
            if not checked(dest(method,gate)):continue
            gate_record=state['jobs'][f'{method}/{gate}']
            if gate_record.get('status')!='complete' or not gate_record.get('accounting'):continue
            measured=[]
            for row in gate_record['accounting']:
                if row.get('MaxRSS'):
                    match=re.fullmatch(r'([0-9.]+)([KMGTP]?)',row['MaxRSS'])
                    assert match,row['MaxRSS']
                    scale={'':1,'K':1024,'M':1024**2,'G':1024**3,'T':1024**4,'P':1024**5}[match[2]]
                    measured.append(float(match[1])*scale/1024**3)
            if not measured:continue
            if slots>=5:continue
            # Explicit quadratic headroom for original R HDBSCAN; actual peak is
            # reported later, never inferred from these reservation amounts.
            lower={30000:128,50000:192,100000:384,120000:480} if method=='DG-scRNA' else {30000:64,50000:96,100000:192,120000:256}
            estimate=max(measured)*(n/gate)**(2 if method=='DG-scRNA' else 1)*1.5
            # nextgen nodes advertise515456MiB. A512GiB reservation cannot fit.
            mem=min(480,max(lower[n],int((estimate+31)//32)*32))
            wall='24:00:00' if method=='DG-scRNA' and n>=100000 else '12:00:00' if method=='DG-scRNA' else '08:00:00'
            script='scalability_job.py' if method=='DG-scRNA' else 'scalability_method_job.py'
            args=[n] if method=='DG-scRNA' else [method,n]
            job=submit(source,script,args,[f'--job-name=claim_scale_{method}_{n}','--cpus-per-task=4',f'--mem={mem}G',f'--time={wall}'])
            record.update(job=job,submitted=time.time(),status='submitted',source=str(source),memory_GB=mem,
                          pilot_gate=gate,gate_maxRSS_GB=max(measured),uncapped_estimate_GB=estimate);slots+=1
            write_json(path,state)
    completed=sum(r.get('status')=='complete' for r in state['jobs'].values())
    state.update(expected=15,completed=completed,last_check=utc(),
                 review_required={k:r for k,r in state['jobs'].items() if r.get('status')=='needs_review'})
    write_json(path,state);print('SCALABILITY_PROGRESS',completed,'/15',flush=True)
    return completed==15

if __name__=='__main__':
    while not tick():time.sleep(45)
