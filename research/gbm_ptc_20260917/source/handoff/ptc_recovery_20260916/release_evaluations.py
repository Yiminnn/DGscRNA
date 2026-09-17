"""Operational completion-gated scheduling; reads no matrices or evaluation labels."""
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
from ptc_common import BASE,GROUPS,geometry_dir,task_list,write_json,utc
assert os.environ.get('SLURM_JOB_ID')
job=os.environ['PTC_EVALUATION_JOB']
assert (BASE/'evaluations/task_000/COMPLETE').exists()
tasks=task_list()
directories=[geometry_dir(t) for t in tasks]
directories += [BASE/'pooled'/g/c for g in GROUPS for c in ['NONE','CCA','HARMONY']]
directories += [BASE/'matched_samples'/s for ss in GROUPS.values() for s in ss]
assert len(directories)==510
eligible=set(json.loads(Path(os.environ['PTC_EVALUATION_TASK_LIST']).read_text()))
released=set();failures=[]
record=BASE/'operational/evaluation_dependency_releases.json'
if record.exists():
 old=json.loads(record.read_text());assert old['evaluation_job']==job;released=set(old['released'])
def release(i):
 r=subprocess.run(['scontrol','release',f'{job}_{i}'],capture_output=True,text=True)
 return i,r.returncode,r.stdout+r.stderr
while not (BASE/'operational/STOP_EVALUATION_CONTROLLER').exists():
 ready=[i for i,d in enumerate(directories) if i in eligible and i not in released and (d/'ANNOTATION_COMPLETE').exists()]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
  for i,code,msg in ex.map(release,ready):
   if code==0:released.add(i)
   else:failures.append(dict(task=i,code=code,message=msg,time=utc()))
 write_json(record,dict(evaluation_job=job,released=sorted(released),failures=failures,
                       updated_at=utc(),controller_job=os.environ['SLURM_JOB_ID']))
 print(utc(),len(released),'evaluation tasks released',flush=True)
 if len(released)==len(eligible):break
 time.sleep(30)
