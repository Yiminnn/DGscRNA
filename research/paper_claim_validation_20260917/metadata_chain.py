"""Short metadata allocations renew between ticks, without stopping science jobs."""
import importlib
import json
import os
from pathlib import Path
import sys
import time
from common import OUT,require_slurm,write_json,utc
from dispatch import freeze,submit

ALLOWED={'analysis_dispatch','full_GBM_launcher','scalability_dispatch','scina_library_dispatch','aux_dispatch','ptc_dispatch'}

def run(module_name):
    require_slurm();assert module_name in ALLOWED
    started=time.monotonic();module=importlib.import_module(module_name)
    path=OUT/'metadata_chains'/(module_name+'.json')
    state=json.loads(path.read_text()) if path.exists() else dict(history=[])
    job=os.environ['SLURM_JOB_ID']
    record=dict(job=job,source=str(Path(__file__).resolve().parent),started_at=utc(),status='running')
    state['history'].append(record);state['active_job']=job;write_json(path,state)
    while True:
        done=module.tick()
        if done:
            record.update(status='completed',finished_at=utc());state['status']='completed'
            write_json(path,state);return
        if time.monotonic()-started>=50*60:
            # Renew only between complete ticks, after all submissions/state writes.
            # These are metadata observers; their submitted fitting jobs are untouched.
            source=freeze('metadata_chain')
            next_job=submit(source,'metadata_chain.py',[module_name],[f'--job-name=claim_meta_{module_name}',
                '--cpus-per-task=1','--mem=2G','--time=01:00:00',f'--dependency=afterany:{job}'])
            record.update(status='renewed_between_ticks',next_job=next_job,finished_at=utc())
            state.update(status='awaiting_next_allocation',next_job=next_job);write_json(path,state);return
        time.sleep(45)

if __name__=='__main__':run(sys.argv[1])
