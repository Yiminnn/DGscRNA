"""Invoke isolated Python3.8 predictor, then apply the frozen L1 evaluation."""
import json
import os
from pathlib import Path
import subprocess
import sys
from common import OUT, require_slurm
if __name__=='__main__':
    require_slurm()
    import evaluate_comparator
    if len(sys.argv)>1 and sys.argv[1]=='cohort':
        import csv
        sample=list(csv.DictReader((OUT/'protocol/sample_order.csv').open()))[int(os.environ['SLURM_ARRAY_TASK_ID'])]['sample']
    elif len(sys.argv)>1:sample=sys.argv[1]
    else:
        sample=(OUT/'protocol/pilot_samples.txt').read_text().split()[int(os.environ['SLURM_ARRAY_TASK_ID'])]
    source=Path(__file__).resolve().parent
    env=os.environ.copy();env['DGLBACKEND']='pytorch'
    subprocess.run([str(OUT/'vendor_envs/scDeepSort/bin/python'),str(source/'deepsort_predict.py'),sample],check=True,env=env)
    evaluate_comparator.run('scDeepSort',sample)
