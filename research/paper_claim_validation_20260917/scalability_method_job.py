"""Cold-input comparison methods; total timing includes required input conversion."""
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
from common import OUT, RSCRIPT, require_slurm, checked, complete, sha, write_json, utc

def run(method,n,repeat=0):
    require_slurm();n=int(n);repeat=int(repeat)
    from scalability_paths import output,input_alias
    sample=input_alias(method,n,repeat);dest=output(method,n,repeat);dest.mkdir(parents=True,exist_ok=True)
    if checked(dest):return
    source=Path(__file__).resolve().parent
    start=time.monotonic()
    if method=='SCINA':
        assert not (dest/'fit_manifest.json').exists(),'Preserve partial attempts before a cold-input repeat'
        subprocess.run([RSCRIPT,str(source/'scalability_SCINA.R'),str(n),sample,str(dest)],check=True)
        fit=json.loads((dest/'fit_manifest.json').read_text());assert fit['status']=='completed'
        model_scope='Frozen glioma marker library; default overlap removal and Unknown enabled'
    elif method=='scDeepSort':
        assert not (OUT/'comparators/scDeepSort'/sample/'COMPLETE').exists(),'Do not time a cached prediction as cold'
        install=json.loads((OUT/'scDeepSort_install.json').read_text())
        assert install['status']=='installed_and_checkpoint_loaded'
        python=str(OUT/'vendor_envs/scDeepSort/bin/python')
        subprocess.run([python,str(source/'deepsort_predict.py'),sample],check=True)
        fit=json.loads((OUT/'comparators/scDeepSort'/sample/'manifest.json').read_text())
        model_scope='Published human Brain atlas model; fixed pretraining information differs from marker methods'
    else:raise ValueError(method)
    elapsed=time.monotonic()-start
    write_json(dest/'manifest.json',dict(status='completed',method=method,n_cells=n,resource_repeat=repeat,input_sample=sample,
        pipeline_elapsed_seconds=elapsed,cold_input=True,input_conversion_included=True,
        scope=model_scope,resource_only=True,accuracy_not_evaluated_on_pooled_input=True,
        hardware=platform.node(),SLURM_CPUS_PER_TASK=os.environ.get('SLURM_CPUS_PER_TASK'),
        fit_manifest=fit,job=os.environ['SLURM_JOB_ID'],source_sha256=sha(__file__),completed_at=utc()))
    complete(dest)

if __name__=='__main__':run(*sys.argv[1:])
