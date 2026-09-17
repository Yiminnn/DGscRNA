"""Operational controller: release only ready tasks after verified global pilot gates."""
from pathlib import Path
import concurrent.futures
import json
import os
import subprocess
import time
from ptc_common import BASE,task_list,geometry_dir,write_json,sha,utc
assert os.environ.get('SLURM_JOB_ID')
annotation_job='7331316'
gate=BASE/'single_sample/MT-1/all__direct__PCA2__s42/HDBSCAN/score_RNA/REFINEMENT_COMPLETE'
assert gate.exists()
assert (BASE/'verification/implementation/COMPLETE').exists()
assert (BASE/'verification/legacy_density/COMPLETE').exists()
tasks=task_list();released=set();failed=[]
path=BASE/'operational/annotation_dependency_releases.json'
if path.exists():released=set(json.loads(path.read_text())['released'])
def release(i):
    result=subprocess.run(['scontrol','update',f'JobId={annotation_job}_{i}','Dependency='],capture_output=True,text=True)
    return i,result.returncode,result.stdout+result.stderr
while not (BASE/'operational/STOP_RELEASE_CONTROLLER').exists():
    ready=[]
    for t in tasks:
        i=t['task_id']
        if i in released:continue
        dest=geometry_dir(t)
        if not (dest/'GEOMETRY_COMPLETE').exists():continue
        assert (dest/'GEOMETRY_COMPLETE').read_text().strip()==sha(dest/'geometry_manifest.json')
        m=json.loads((dest/'geometry_manifest.json').read_text())
        assert m['status']=='completed' and m['task']==t
        ready.append(i)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        for i,code,msg in ex.map(release,ready):
            if code==0:released.add(i)
            else:failed.append(dict(task=i,code=code,message=msg,time=utc()))
    write_json(path,dict(annotation_job=annotation_job,released=sorted(released),failures=failed,
                        updated_at=utc(),controller_job=os.environ['SLURM_JOB_ID']))
    print(utc(),len(released),'dependencies released',flush=True)
    if len(released)==len(tasks):break
    time.sleep(30)
