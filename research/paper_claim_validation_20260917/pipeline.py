"""One SLURM allocation, sequential native R stages and explicit terminal DL."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from common import ROOT, OUT, RSCRIPT, PYTHON, FEATURES, ROUTES, require_slurm, write_json, sha, utc

def run_stage(argv,env=None):
    print('EXEC',json.dumps(argv),flush=True)
    subprocess.run(argv,check=True,env=env,cwd=ROOT)

def unit(sample,budget,routes,pilot=False,seed=42):
    require_slurm()
    started=time.monotonic()
    source=Path(__file__).resolve().parent
    condition=budget if seed==42 else f'{budget}_seed{seed}'
    prep=OUT/'GBM'/sample/condition
    prep.mkdir(parents=True,exist_ok=True)
    run_stage([RSCRIPT,str(source/'prepare_R.R'),sample,budget,str(seed)])
    for route in routes:
        env=os.environ.copy();env['DGSCRNA_ONLY_ROUTE']=route
        run_stage([RSCRIPT,str(source/'score_R.R'),sample,str(prep)],env)
        run_stage([PYTHON,'-u',str(source/'terminal.py'),str(prep/route),'L00_mean' if pilot else 'all'])
    receipt=dict(status='pilot_terminal_complete' if pilot else 'unit_terminal_complete',sample=sample,budget=budget,
        routes=routes,seed=seed,source_bundle=str(source),driver_sha256=sha(__file__),
        job=os.environ['SLURM_JOB_ID'],elapsed_seconds=time.monotonic()-started,completed_at=utc())
    write_json(prep/('pilot_manifest.json' if pilot else 'pipeline_manifest.json'),receipt)
    print('PIPELINE_COMPLETE',sample,budget,pilot,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['pilot','unit','tasklist'])
    p.add_argument('rest',nargs='*');a=p.parse_args()
    if a.mode=='pilot':
        samples=(OUT/'protocol/pilot_samples.txt').read_text().split()
        unit(samples[int(os.environ['SLURM_ARRAY_TASK_ID'])],'hvg2000',['PCA30_SNN','UMAP2_HDBSCAN_R'],pilot=True)
    elif a.mode=='unit':unit(a.rest[0],a.rest[1],ROUTES)
    else:
        t=json.loads(Path(a.rest[0]).read_text())[int(os.environ['SLURM_ARRAY_TASK_ID'])]
        unit(t['sample'],t['budget'],t.get('routes',ROUTES),seed=t.get('seed',42))
