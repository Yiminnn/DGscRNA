"""Administrative A2 status from explicit job handles and completion manifests."""
from datetime import datetime,timezone
import json
from pathlib import Path
import subprocess

ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/no_clustering'
registry=json.loads((OUT/'jobs.json').read_text()) if (OUT/'jobs.json').exists() else []
for job in registry:
    query=subprocess.run(['sacct','-n','-X','-P','-j',str(job['job_id']),'--format=JobIDRaw,State,ExitCode'],text=True,capture_output=True,timeout=15)
    job['observed_jobs']=[line.split('|')[:3] for line in query.stdout.splitlines() if line.strip()]
    job['observation_status']='ok' if query.returncode==0 else 'unknown'
units=[]
for path in (OUT/'GBM').glob('*/*/evaluation/COMPLETE'):
    units.append('/'.join(path.parts[-4:-2]))
summary_ready=(OUT/'summary/COMPLETE').exists()
status=dict(stage='A2_GBM',work_package='A',status='completed_pending_parent_validation' if summary_ready else 'running',
    updated_at=datetime.now(timezone.utc).isoformat(),jobs=registry,
    summary='GBM cellwise seed-mechanism replacement; original marker and DL inputs preserved. No PTC computation.',
    completed_evaluation_units=len(units),expected_evaluation_units=242,
    expected_candidate_conditions=1210,completed_units=units,
    original_regression_manifest_present=(OUT/'regression/COMPLETE').exists(),
    aggregate_manifest_present=summary_ready,whole_work_package_A_complete=False,
    evidence=[str((OUT/name).relative_to(ROOT)) for name in ['protocol.json','source/PROTOCOL.md','source/SOURCE_MANIFEST.json','jobs.json']],
    limitations=['Explicit seed construction replacement, not pure deletion of clustering','All five lambda values retained','Patient-based retrospective label-heldout selection','No-op and unavailable trained refinement reported explicitly'])
path=OUT/'status.json';temp=path.with_suffix('.tmp');temp.write_text(json.dumps(status,ensure_ascii=False,indent=2)+'\n');temp.replace(path)
print(json.dumps(dict(status=status['status'],completed_evaluation_units=len(units),jobs=len(registry))))
