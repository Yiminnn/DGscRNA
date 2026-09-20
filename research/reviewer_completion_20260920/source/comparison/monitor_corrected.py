"""Read-only SLURM/status monitor. Stops on failed task or final summary terminal state."""
import json,subprocess,time,sys
from pathlib import Path
root=Path('/fs/scratch/PCON0080/yimin/dgscrna')
status=root/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/comparison/scDeepSort_LogNormalize/status.json'
last=None
while True:
    account=subprocess.run(['sacct','-j','7446013,7446014','--noheader','--parsable2','--format=JobID,State,ExitCode'],text=True,capture_output=True)
    if account.returncode:
        print('TRANSIENT_SACCT',account.stderr.strip(),flush=True);time.sleep(30);continue
    rows=[x.split('|') for x in account.stdout.splitlines() if x.strip()]
    roots=[r for r in rows if '.' not in r[0]]
    bad=[r for r in roots if any(x in r[1] for x in ['FAILED','OUT_OF_MEMORY','TIMEOUT','CANCELLED','NODE_FAIL','BOOT_FAIL'])]
    if bad:
        print('ACTION_REQUIRED',json.dumps(bad),flush=True);sys.exit(2)
    state=json.loads(status.read_text()) if status.exists() else {}
    count=state.get('n_samples_complete',3)
    summary=[r for r in roots if r[0]=='7446014']
    stamp=(count,summary[0][1] if summary else 'unobserved')
    if stamp!=last:
        print('COHORT_PROGRESS',json.dumps(dict(n_complete=count,total=121,summary=summary)),flush=True);last=stamp
    if summary and summary[0][1]=='COMPLETED':
        print('SUMMARY_TERMINAL',json.dumps(summary),flush=True);break
    # Pending/running job identity remains authoritative; no automatic restart.
    time.sleep(30)
