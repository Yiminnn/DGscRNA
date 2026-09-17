#!/usr/bin/env python3
"""Submit dependent GBM delivery jobs after all input gates close; no fitting here."""
import json
import subprocess
import time
from common import ROOT,OUT,samples,utc,write_json,require_slurm
require_slurm()
ledger=OUT/'finalization_ledger.json'
state=json.loads(ledger.read_text()) if ledger.exists() else {'jobs':{}}
scripts=ROOT/'handoff/hvg_ptc_20260916'
stage_rules={
    'summary':('summarize.sbatch',lambda:all((OUT/'evaluation'/s/'COMPLETE').exists() and (OUT/'quality'/s/'COMPLETE').exists() for s in samples()),OUT/'summary/COMPLETE'),
    'verify':('verify_results.sbatch',lambda:all((OUT/'evaluation'/s/'COMPLETE').exists() for s in samples()),OUT/'verification/COMPLETE'),
    'secondary_verify':('verify_singleton.sbatch',lambda:all((OUT/'singleton_robustness'/s/'COMPLETE').exists() for s in samples()),OUT/'singleton_robustness_summary/COMPLETE'),
    'methods':('compare_methods.sbatch',lambda:done('summary') and done('secondary_verify'),OUT/'method_comparisons/COMPLETE'),
    'plots':('plots.sbatch',lambda:done('methods'),OUT/'figures/figure_manifest.json'),
    'notebook':('build_notebook.sbatch',lambda:done('plots') and done('verify'),OUT/'notebooks/execution_manifest.json'),
    'report':('report.sbatch',lambda:done('notebook'),OUT/'GBM_DELIVERY.json'),
}
def done(stage):return state['jobs'].get(stage,{}).get('state')=='COMPLETED'
while True:
    for stage,(script,ready,marker) in stage_rules.items():
        if done(stage):continue
        row=state['jobs'].get(stage)
        if row:
            p=subprocess.run(['sacct','-X','-n','-P','-j',row['job'],'--format=JobIDRaw,State,ExitCode'],capture_output=True,text=True,check=True)
            records=[line.split('|') for line in p.stdout.splitlines() if line.split('|')[0]==row['job']]
            if records:
                status=records[0][1];row['slurm_state']=status
                if status=='COMPLETED':
                    assert marker.exists(),(stage,'Missing delivery artifact',str(marker))
                    row.update(state='COMPLETED',completed_at=utc())
                    print(utc(),'completed',stage,row['job'],flush=True)
                elif status not in ['RUNNING','PENDING','COMPLETING','CONFIGURING']:
                    row.update(state='FAILED',failed_at=utc())
                    write_json(ledger,state)
                    raise RuntimeError(f'{stage} job {row["job"]}: {status}; inspect and retry explicitly')
            continue
        if ready():
            job=subprocess.run(['sbatch','--parsable',str(scripts/script)],capture_output=True,text=True,check=True).stdout.strip().split(';')[0]
            state['jobs'][stage]=dict(job=job,state='SUBMITTED',submitted_at=utc())
            print(utc(),'submitted',stage,job,flush=True)
            write_json(ledger,state)
    state['updated_at']=utc();write_json(ledger,state)
    if done('report'):break
    time.sleep(30)
