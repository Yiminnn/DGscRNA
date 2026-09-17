import datetime,hashlib,json,os,shutil
from pathlib import Path
from stage_delivery import ROOT,OUT,PARENT,STAGE,REMOTE,sha

assert os.environ.get('SLURM_JOB_ID')
plan=PARENT/'R_reference_campaign_20260917_delivery_manifest.json'
m=json.loads(plan.read_text())
record=dict(status='delivered_and_verified',time=datetime.datetime.now(datetime.timezone.utc).isoformat(),
  job=os.environ['SLURM_JOB_ID'],remote=REMOTE,n_files=m['n_files'],total_bytes=m['total_bytes'],
  manifest_sha256=sha(plan),notebook_sha256=sha(ROOT/'notebooks/dgscrna_results.ipynb'),
  verification='rclone copy and one-way full-download check both exited 0; no differences',
  old_results_preserved=True,website_modified=False,new_version_notebook_created=False)
p=OUT/'summary/DELIVERY_RECEIPT.json';p.write_text(json.dumps(record,indent=2)+'\n')
target=STAGE/p.relative_to(ROOT);shutil.copy2(p,target)
print(json.dumps(record,indent=2),flush=True)
