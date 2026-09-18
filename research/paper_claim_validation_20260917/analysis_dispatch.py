"""Start independent GBM summaries as their own prerequisites finish."""
import csv
import json
import subprocess
import time
from common import OUT,checked,write_json,utc
from dispatch import freeze,submit,states

def tick():
    path=OUT/'analysis_dispatch_state.json'
    state=json.loads(path.read_text()) if path.exists() else dict(jobs={})
    active=states();aux=json.loads((OUT/'aux_dispatch_state.json').read_text())
    rep=json.loads((OUT/'representation_dispatch_state.json').read_text())
    scale=json.loads((OUT/'scalability_dispatch_state.json').read_text())
    cohort=list(csv.DictReader((OUT/'protocol/cohort.csv').open()))
    core=checked(OUT/'summary','aggregate_manifest.json','AGGREGATE_COMPLETE')
    controls=all(sum(k.startswith(prefix) and v.get('status')=='complete' for k,v in aux['jobs'].items())==n
        for prefix,n in [('geometry/',363),('DL_control/',54)]) and rep.get('completed')==180
    comparators=all(checked(OUT/'comparators'/m/r['sample']/'evaluation')
        for m in ['scType','scCATCH','SCINA','SingleR','scDeepSort'] for r in cohort)
    compared=checked(OUT/'comparison_summary')
    selected=checked(OUT/'DG_fixed_partition_selection')
    tasks=[('DG_fixed_partition_selection','fixed_partition_selection.py',core,'8G','00:30:00'),
        ('controls_summary','aggregate_controls.py',core and controls,'32G','02:00:00'),
        ('workflow_choice_summary','workflow_choice_summary.py',core and checked(OUT/'controls_summary'),'8G','00:30:00'),
        ('comparison_summary','aggregate_comparators.py',core and comparators and selected,'32G','02:00:00'),
        ('comparison_summary/diagnostics','comparison_diagnostics.py',compared,'24G','02:00:00'),
        ('unknown_summary','unknown_diagnostics.py',selected,'16G','02:00:00'),
        ('unknown_expression','unknown_expression.py',selected,'32G','04:00:00'),
        ('scalability_summary','aggregate_scalability.py',scale.get('completed')==45,'8G','01:00:00')]
    slots=sum(name.startswith('claim_GBM_analysis_') for name,status in active.values())
    for folder,script,ready,mem,wall in tasks:
        r=state['jobs'].setdefault(folder,{})
        if checked(OUT/folder):r['status']='complete';continue
        if not ready or r.get('job') in active:continue
        if r.get('job') and time.time()-r['submitted']<180:continue
        if r.get('job'):
            raw=subprocess.check_output(['sacct','-j',r['job'],'-X','-n','-P','--format=JobID,State,ExitCode'],text=True)
            done=[s.split('|') for s in raw.splitlines() if s.startswith(r['job']+'|')]
            if done and done[0][1].split()[0] in ['FAILED','TIMEOUT','OUT_OF_MEMORY','COMPLETED','CANCELLED','NODE_FAIL']:
                r.update(status='needs_review',accounting=done[0])
            continue
        if slots>=3:continue
        source=freeze('GBM_analysis')
        r.update(job=submit(source,script,[],[f'--job-name=claim_GBM_analysis_{script[:-3]}',
            '--cpus-per-task=4',f'--mem={mem}',f'--time={wall}']),source=str(source),submitted=time.time(),status='submitted')
        slots+=1;write_json(path,state)
    state.update(expected=len(tasks),completed=sum(r.get('status')=='complete' for r in state['jobs'].values()),last_check=utc(),
        review_required={k:v for k,v in state['jobs'].items() if v.get('status')=='needs_review'})
    write_json(path,state);print('ANALYSIS_PROGRESS',state['completed'],'/',len(tasks),flush=True)
    return state['completed']==len(tasks)

if __name__=='__main__':
    while not tick():time.sleep(45)
