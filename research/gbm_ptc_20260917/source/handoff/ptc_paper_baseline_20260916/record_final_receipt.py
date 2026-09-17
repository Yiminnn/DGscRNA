from pathlib import Path
import os,json,hashlib
from datetime import datetime,timezone
p=Path('/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916')
m=json.loads((p/'PTC_baseline_correction_staging.json').read_text())
log=p/'PTC_baseline_correction_final_remote_check.log'
s=log.read_text()
assert '0 differences found' in s and f"{m['n_files']} matching files" in s
m.update(completed_utc=datetime.now(timezone.utc).isoformat(),job=os.environ.get('SLURM_JOB_ID'),
    verification='rclone copy and full download check exit 0; all listed files match; zero differences',
    check_log_sha256=hashlib.sha256(log.read_bytes()).hexdigest())
(p/'PTC_baseline_correction_upload_receipt_final.json').write_text(json.dumps(m,indent=2)+'\n')
print(f"Final correction delivery verified: {m['n_files']} files, zero differences.")
