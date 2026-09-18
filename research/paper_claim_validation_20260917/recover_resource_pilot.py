"""One bounded pilot retry for a confirmed resource failure, preserving attempts."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from common import OUT,RSCRIPT,require_slurm,checked,sha,write_json,utc

def run(method,previous_job):
    require_slurm()
    if method=='SCINA':
        dest=OUT/'comparators/SCINA/SN040/L15';manifest='manifest.json';flag='COMPLETE'
        script='scina_R.R';args=['SN040','15']
    elif method=='scCATCH':
        dest=OUT/'comparators/scCATCH/SN040/hvg2000/UMAP2_HDBSCAN_R'
        manifest='cohort_manifest.json';flag='COHORT_COMPLETE';script='sccatch_R.R'
        args=['SN040','cohort','hvg2000','UMAP2_HDBSCAN_R']
    else:raise ValueError(method)
    if checked(dest,manifest,flag):
        print('PILOT_ALREADY_VERIFIED_NO_RETRY',method,flush=True);return
    assert not (dest/flag).exists(),'An invalid completion flag requires manual review'
    records=[]
    for attempt in range(15):
        raw=subprocess.check_output(['sacct','-j',previous_job,'-X','-n','-P','--format=JobID,State,ExitCode'],text=True)
        records=[r.split('|') for r in raw.splitlines() if r.startswith(previous_job+'|')]
        if records and records[0][1].split()[0] in ['TIMEOUT','OUT_OF_MEMORY','NODE_FAIL']:break
        if records and records[0][1].split()[0] in ['FAILED','CANCELLED','COMPLETED']:
            raise RuntimeError('Not a verified resource failure; inspect the original attempt: '+repr(records))
        time.sleep(4)
    assert records and records[0][1].split()[0] in ['TIMEOUT','OUT_OF_MEMORY','NODE_FAIL'],records
    archive=dest/'attempts'/('resource_failure_'+previous_job);archive.mkdir(parents=True,exist_ok=True)
    files={}
    for p in dest.iterdir():
        if not p.is_file():continue
        target=archive/p.name
        if target.exists():assert sha(target)==sha(p)
        else:shutil.copy2(p,target)
        files[p.name]=sha(target)
    write_json(archive/'RESOURCE_RECOVERY.json',dict(previous_job=previous_job,accounting=records,
        retry_job=os.environ['SLURM_JOB_ID'],scientific_parameters_unchanged=True,preserved_files=files,
        decision='One larger-walltime allocation only after a scheduler-confirmed resource failure',time=utc()))
    source=Path(__file__).resolve().parent
    subprocess.run([RSCRIPT,str(source/script),*args],check=True)
    assert checked(dest,manifest,flag)

if __name__=='__main__':run(*sys.argv[1:])
