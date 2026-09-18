"""Metadata-only resumable SLURM dispatcher. Never evaluates matrices locally."""
import argparse
import csv
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
from common import ROOT, CODE, OUT, FEATURES, ROUTES, sha, write_json, checked, utc

def freeze(label):
    files=[p for p in CODE.iterdir() if p.suffix in ['.py','.R','.sbatch']]
    import hashlib
    sig=hashlib.sha256(''.join(p.name+sha(p) for p in sorted(files)).encode()).hexdigest()[:16]
    dest=OUT/'source_snapshots'/f'{label}_{sig}'
    dest.mkdir(parents=True,exist_ok=True)
    for p in files:
        if (dest/p.name).exists():assert sha(dest/p.name)==sha(p)
        else:shutil.copy2(p,dest/p.name)
    write_json(dest/'SOURCE_MANIFEST.json',{p.name:sha(p) for p in files})
    return dest

def submit(source,script,args=(),options=()):
    cmd=['sbatch','--parsable',*options,str(source/'job.sbatch'),str(source/script),*map(str,args)]
    p=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True,check=True)
    jid=p.stdout.strip().split(';')[0]
    with (OUT/'submissions.jsonl').open('a') as f:f.write(json.dumps(dict(time=utc(),job=jid,command=cmd))+'\n')
    print('SUBMITTED',jid,script,flush=True)
    return jid

def states():
    lines=subprocess.check_output(['squeue','-u','yimin','-h','-r','-o','%i|%j|%T'],text=True).splitlines()
    return {r.split('|')[0]:r.split('|')[1:] for r in lines}

def tick():
    OUT.mkdir(parents=True,exist_ok=True)
    with (OUT/'dispatch.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        path=OUT/'dispatch_state.json'
        state=json.loads(path.read_text()) if path.exists() else dict(phase='waiting_for_verified_pilots',units={},attempts={})
        if not checked(OUT/'verification','pilot_manifest.json','PILOT_PASSED'):
            state['last_check']=utc();write_json(path,state);return False
        for name in ['SN040_cluster_rng_audit.json','SN040_wilcox_kernel_audit.json']:
            p=OUT/'verification'/name
            if not p.exists() or json.loads(p.read_text()).get('status')!='passed':
                state['phase']='waiting_for_pilot_rng_checks';state['last_check']=utc();write_json(path,state);return False
        active=states();source=freeze('core')
        samples=list(csv.DictReader((OUT/'protocol/sample_order.csv').open()))
        ordered=['hvg2000','hvg5000','hvg500','hvg1000','hvg3000','all']
        tasks=[dict(sample=s['sample'],budget=b) for s in samples for b in ordered]
        taskfile=OUT/'protocol/core_tasks.json'
        if not taskfile.exists():write_json(taskfile,tasks)
        else:assert json.loads(taskfile.read_text())==tasks
        if not state.get('main_array'):
            state['main_array']=submit(source,'pipeline.py',['tasklist',taskfile],
                ['--job-name=claim_GBM_core',f'--array=0-{len(tasks)-1}%64','--cpus-per-task=4','--mem=40G','--time=04:00:00'])
            state['main_source']=str(source);state['main_submitted']=time.time();state['phase']='GBM_core_running';write_json(path,state)
            # Account for the just-submitted tasks immediately; scheduler queries may lag.
            active=states()
        source=Path(state.get('post_source',state['main_source']))
        active_post=sum(v[0]=='claim_GBM_post' for v in active.values())
        active_retry=sum(v[0]=='claim_GBM_retry' for v in active.values())
        n_complete=0;n_terminal=0
        for i,t in enumerate(tasks):
            key=t['sample']+'/'+t['budget'];u=state['units'].setdefault(key,{})
            prep=OUT/'GBM'/t['sample']/t['budget']
            fit_complete=(prep/'pipeline_manifest.json').exists()
            if fit_complete:
                n_terminal+=1
                if checked(prep/'evaluation') and checked(prep/'figures','manifest.json','FIGURES_COMPLETE'):
                    u['status']='complete';n_complete+=1;continue
                if not checked(OUT/'GBM'/t['sample']/'hvg2000','prepare_manifest.json','PREPARED'):continue
                jid=u.get('post_job')
                if jid and jid in active:continue
                if jid and jid not in active and time.time()-u.get('post_submitted',0)<180:continue
                tries=u.get('post_attempts',0)
                if tries>=3:
                    u['needs_review']='postprocess_failed';continue
                if active_post<24:
                    u['post_job']=submit(source,'postprocess.py',[prep],['--job-name=claim_GBM_post','--cpus-per-task=2','--mem=12G','--time=00:30:00'])
                    u['post_attempts']=tries+1;u['post_submitted']=time.time();active_post+=1;write_json(path,state)
                continue
            job=f'{state["main_array"]}_{i}'
            if job in active:continue
            retry=u.get('retry_job')
            if retry and retry in active:continue
            if not retry and time.time()-state.get('main_submitted',0)<180:continue
            # Never infer failure from a brief scheduler gap; wait for an actual terminal accounting state.
            jid=retry or job
            raw=subprocess.run(['sacct','-j',jid,'-X','-n','-P','--format=JobID,State,ExitCode'],text=True,capture_output=True,check=True).stdout
            lines=[r.split('|') for r in raw.splitlines() if r.strip()]
            terminal=[r for r in lines if r[0]==jid and r[1].split()[0] in ['FAILED','TIMEOUT','OUT_OF_MEMORY','CANCELLED','NODE_FAIL','PREEMPTED','COMPLETED']]
            if not terminal:continue
            u.setdefault('failures',[])
            if not any(r['job']==jid for r in u['failures']):u['failures'].append(dict(job=jid,state=terminal[0][1],exit=terminal[0][2]))
            tries=u.get('retry_attempts',0)
            # Retry OOM/time/node failures with more resources; algorithm errors require inspection.
            why=terminal[0][1]
            if tries>=2 or not any(s in why for s in ['TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED']):
                u['needs_review']='fit_failed:'+why;continue
            if active_retry<4:
                u['retry_job']=submit(source,'pipeline.py',['unit',t['sample'],t['budget']],
                    ['--job-name=claim_GBM_retry','--cpus-per-task=4','--mem=80G','--time=08:00:00'])
                u['retry_attempts']=tries+1;active_retry+=1;write_json(path,state)
        state.update(last_check=utc(),n_expected=len(tasks),n_terminal_complete=n_terminal,n_evaluated_plotted=n_complete,
            review_required={k:u['needs_review'] for k,u in state['units'].items() if 'needs_review' in u})
        if n_complete==len(tasks):state['phase']='GBM_core_complete_followups_pending'
        write_json(path,state)
        print('PROGRESS',n_terminal,n_complete,'of',len(tasks),'review',len(state['review_required']),flush=True)
        return n_complete==len(tasks)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--loop',action='store_true');a=p.parse_args()
    while True:
        done=tick()
        if done or not a.loop:break
        time.sleep(45)
