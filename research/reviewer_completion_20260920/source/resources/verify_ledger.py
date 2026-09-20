from pathlib import Path
from collections import Counter
import csv,json,math,os,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from collect_resources import OUT,sha,write_json,utc,seconds
assert os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_STEP_ID')
snapshot=Path(json.loads((OUT/'latest.json').read_text())['snapshot']);m=json.loads((snapshot/'manifest.json').read_text())
for name,digest in m['files'].items():assert sha(snapshot/name)==digest,(name,'hash')
assert json.loads((snapshot/'missing_accounting_job_ids.json').read_text())==[]
raw=list(csv.DictReader((snapshot/'accounting_raw.csv').open()));rawmap={(r['DBIndex'],r['JobID']):r for r in raw}
rows=list(csv.DictReader((snapshot/'job_attempt_ledger.csv').open()));assert len({(r['db_index'],r['job_id']) for r in rows})==len(rows)
keys={(r['db_index'],r['job_id']) for r in rows if r['shared_allocation_excluded_from_totals']=='False'}
for r in rows:
 assert r['user']=='yimin' and r['job_id']!='7204040'
 a=rawmap[(r['db_index'],r['job_id'])]
 assert math.isclose(float(r['allocated_cpu_hours']),float(a['ElapsedRaw'])*int(a['AllocCPUS'])/3600,abs_tol=1e-10)
 assert math.isclose(float(r['allocated_cpu_hours']),float(a['CPUTimeRAW'])/3600,abs_tol=1e-10)
 if r['measured_cpu_hours']:assert math.isclose(float(r['measured_cpu_hours']),seconds(a['TotalCPU'])/3600,abs_tol=1e-10)
 if '.' in r['job_id']:assert (r['db_index'],r['job_id'].split('.')[0]) not in keys,'Parent/step double count'
 for f in ['batch_maxrss_kib','max_step_maxrss_kib']:
  if r[f]:assert float(r[f])>=0
assert sum(r['state']=='REQUEUED' for r in rows)>=3
counts=Counter(r['phase'] for r in rows);families=Counter(r['family'] for r in rows)
proof=dict(status='passed_snapshot_ledger_arithmetic_and_coverage',snapshot=str(snapshot),n_rows=len(rows),unique_attempt_keys=True,allocated_CPUTimeRAW_exact=True,measured_TotalCPU_exact=True,no_parent_step_double_count=True,missing_recorded_job_ids=0,identity='yimin',source_sha256=sha(__file__),snapshot_manifest_sha256=sha(snapshot/'manifest.json'),job_step=os.environ['SLURM_JOB_ID']+'.'+os.environ['SLURM_STEP_ID'],completed_at=utc(),phase_counts=dict(counts),family_counts=dict(families),F_closed=False)
write_json(snapshot/'verification.json',proof)
print(json.dumps(proof),flush=True)
