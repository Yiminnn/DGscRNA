"""Release larger cold-input resource runs only after measured smaller pilots."""
import json
import math
import re
import shutil
import subprocess
import time
from common import OUT, checked, write_json, utc
from dispatch import freeze,submit,states
from scalability_paths import output

SIZES=[10000,30000,50000,100000,120000]
REPEATS=[0,1,2]
PILOTS={'DG-scRNA':'7363515','SCINA':'7363795','scDeepSort':'7363796'}

def dest(method,n,repeat=0):return output(method,n,repeat)

def task_key(method,n,repeat=0):return f'{method}/{n}' if repeat==0 else f'{method}/{n}/repeat{repeat}'

def accounting(job):
    lines=subprocess.check_output(['sacct','-j',str(job),'-P','--format=JobID,State,ElapsedRaw,MaxRSS,ExitCode,ReqMem,AllocCPUS,NodeList'],text=True).splitlines()
    data=[dict(zip(lines[0].split('|'),r.split('|'))) for r in lines[1:]]
    return data

def completed_accounting(record):
    rows=record.get('accounting',[])
    parent=[r for r in rows if r.get('JobID')==record.get('job')]
    return (len(parent)==1 and parent[0].get('State')=='COMPLETED'
            and any(r.get('MaxRSS') for r in rows))

def tick():
    path=OUT/'scalability_dispatch_state.json'
    state=json.loads(path.read_text()) if path.exists() else dict(jobs={})
    active=states();source=freeze('scalability')
    slots=sum((name.startswith('claim_scale_') or name.startswith('claim_scalability_')) and name!='claim_scale_dispatch' for name,status in active.values())
    core=json.loads((OUT/'dispatch_state.json').read_text())
    limit=16 if core.get('n_terminal_complete')==726 else 8
    tasks=[(method,n,0) for method in PILOTS for n in SIZES]
    tasks += [(method,n,rep) for n in SIZES for rep in [1,2] for method in ['SCINA','scDeepSort','DG-scRNA']]
    for method,n,repeat in tasks:
            key=task_key(method,n,repeat);record=state['jobs'].setdefault(key,{})
            if n==10000 and repeat==0 and not record:record.update(job=PILOTS[method],submitted=0,status='submitted',pilot=True)
            if checked(dest(method,n,repeat)):
                # A result marker may precede final SLURM accounting. Refresh a
                # partial snapshot instead of persisting RUNNING/empty RSS forever.
                if not completed_accounting(record):record['accounting']=accounting(record['job'])
                exact=[r for r in record['accounting'] if r.get('JobID')==record['job']]
                if completed_accounting(record):record['status']='complete'
                elif exact and exact[0]['State'].split()[0] in ['FAILED','TIMEOUT','OUT_OF_MEMORY','CANCELLED','NODE_FAIL']:
                    record['status']='needs_review'
                else:record['status']='results_ready_awaiting_accounting'
                continue
            if record.get('job') in active:continue
            if record.get('job') and time.time()-record.get('submitted',0)<180:continue
            if record.get('job'):
                acc=accounting(record['job'])
                exact=[r for r in acc if r['JobID']==record['job']]
                if exact and exact[0]['State'].split()[0] in ['FAILED','TIMEOUT','OUT_OF_MEMORY','CANCELLED','COMPLETED','NODE_FAIL']:
                    record.update(status='needs_review',accounting=acc)
                continue
            # The passed10k native-method pilot already gates every <=50k cold
            # input. Independent same-input repetitions need not serialize on
            # the first30k/50k endpoint; preserve its conservative reservation.
            # All100k/120k runs still require the completed measured50k pilot.
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
            if slots>=limit:continue
            # Explicit quadratic headroom for original R HDBSCAN; actual peak is
            # reported later, never inferred from these reservation amounts.
            lower=({10000:64,30000:128,50000:192,100000:384,120000:448} if method=='DG-scRNA' else
                   {10000:64,30000:32,50000:48,100000:64,120000:80} if method=='SCINA' else
                   {10000:64,30000:64,50000:96,100000:192,120000:256})
            estimate=max(measured)*(n/gate)*1.5
            if method=='DG-scRNA':
                # Three additional condensed-distance buffers are reserved above
                # linear scaling of the entire observed pilot footprint. Scaling
                # the full Seurat/DEG heap quadratically over-reserves memory.
                estimate+=3*8*(n*(n-1)-gate*(gate-1))/2/1024**3
            # nextgen also limits usable CPUs to120 and memory/CPU to4027MiB.
            # 480GiB would force123CPUs despite physical RAM fitting;448GiB needs114,
            # leaving room for cores reserved for the installed GPU GRES.
            mem=min(448,max(lower[n],int((estimate+31)//32)*32))
            if repeat:
                same_size=state['jobs'][task_key(method,n)]
                if same_size.get('memory_GB'):mem=same_size['memory_GB']
            wall='12:00:00' if method=='DG-scRNA' else '00:30:00' if method=='SCINA' else '08:00:00'
            wall_details={}
            if method=='DG-scRNA' and n>=100000:
                parent=[r for r in gate_record['accounting'] if r['JobID']==gate_record['job']]
                assert len(parent)==1 and parent[0]['State']=='COMPLETED',parent
                gate_seconds=int(parent[0]['ElapsedRaw']);assert gate_seconds>0
                # Scheduling headroom only, not an empirical complexity model or
                # a reported runtime. Keep every scientific parameter unchanged.
                wall_estimate=gate_seconds*(n/gate)**2*1.5
                hours=min(72,max(24,6*math.ceil(wall_estimate/(6*3600))))
                wall=f'{hours}:00:00'
                wall_details=dict(gate_elapsed_seconds=gate_seconds,
                    uncapped_wall_estimate_seconds=wall_estimate,
                    wall_reservation_policy='1.5 * completed50k_elapsed * (n/50000)^2; round up6h; min24h/max72h')
            script='scalability_job.py' if method=='DG-scRNA' else 'scalability_method_job.py'
            args=[n,repeat] if method=='DG-scRNA' else [method,n,repeat]
            job=submit(source,script,args,[f'--job-name=claim_scale_{method}_{n}_r{repeat}','--cpus-per-task=4',f'--mem={mem}G',f'--time={wall}'])
            record.update(job=job,submitted=time.time(),status='submitted',source=str(source),memory_GB=mem,
                          resource_repeat=repeat,pilot_gate=gate,gate_maxRSS_GB=max(measured),uncapped_estimate_GB=estimate,
                          requested_walltime=wall,**wall_details);slots+=1
            write_json(path,state)
    completed=sum(r.get('status')=='complete' for r in state['jobs'].values())
    state.update(expected=45,repeats_per_method_size=3,completed=completed,last_check=utc(),concurrency_limit=limit,
                 review_required={k:r for k,r in state['jobs'].items() if r.get('status')=='needs_review'})
    write_json(path,state);print('SCALABILITY_PROGRESS',completed,'/45',flush=True)
    return completed==45

if __name__=='__main__':
    while not tick():time.sleep(45)
