"""One sample per SLURM job; cohort progress from checksum-verified completions."""
import csv,fcntl,json,os,sys
from pathlib import Path
from scdeepsort_corrected import run_sample,DEST,ROOT
from common import checked,write_json,utc
if len(sys.argv)>1:
    sample=sys.argv[1]
else:
    roster=(Path(__file__).parent/'scdeepsort_remaining_samples.txt').read_text().split()
    sample=roster[int(os.environ['SLURM_ARRAY_TASK_ID'])]
result=run_sample(sample)
cohort=list(csv.DictReader((ROOT/'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/protocol/cohort.csv').open()))
with (DEST/'progress.lock').open('a+') as lock:
    fcntl.flock(lock,fcntl.LOCK_EX)
    done=[r['sample'] for r in cohort if checked(DEST/r['sample'])]
    task_jobs=sorted({json.loads((DEST/s/'manifest.json').read_text())['job'] for s in done})
    jobs=['7445743','7445842','7445852','7445943','7446013','7446014']
    write_json(DEST/'status.json',dict(stage='C_GBM_SCDEEPSORT',status='predictions_complete_verification_pending' if len(done)==121 else 'running',
        updated_at=utc(),jobs=jobs,array_job='7446013',summary_job='7446014',array_concurrency=16,sample_task_jobs=task_jobs,
        completed=[f'{len(done)}/121 corrected-input predictions; no retraining'],
        remaining=['Complete remaining sample predictions','Full corrected comparison/consensus version and accounting verification'],
        completed_samples=done,n_samples_complete=len(done),n_samples_total=121,
        evidence=[str(DEST/s/'manifest.json') for s in done]))
print('ONE_SAMPLE_COMPLETE',sample,flush=True)
