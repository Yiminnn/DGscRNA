"""Administrative A2 status from explicit job handles and completion manifests."""
from datetime import datetime,timezone
import json
import hashlib
from pathlib import Path
import subprocess

ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/no_clustering'
registry=json.loads((OUT/'jobs.json').read_text()) if (OUT/'jobs.json').exists() else []
for job in registry:
    query=subprocess.run(['sacct','-n','-X','-P','-j',str(job['job_id']),'--format=JobIDRaw,State,ExitCode'],text=True,capture_output=True,timeout=15)
    job['observed_jobs']=[line.split('|')[:3] for line in query.stdout.splitlines() if line.strip()]
    job['observation_status']='ok' if query.returncode==0 else 'unknown'
def manifest_complete(directory):
    try:
        marker=directory/'COMPLETE';manifest=directory/'manifest.json'
        return marker.read_text().strip()==hashlib.sha256(manifest.read_bytes()).hexdigest()
    except OSError:
        return False

units=[];invalid_markers=[]
for path in (OUT/'GBM').glob('*/*/evaluation/COMPLETE'):
    identity='/'.join(path.parts[-4:-2])
    (units if manifest_complete(path.parent) else invalid_markers).append(identity)
units.sort();invalid_markers.sort()
summary_ready=manifest_complete(OUT/'summary')
protocol=json.loads((OUT/'protocol.json').read_text())
active_source=next(path for path in OUT.glob('source*') if path.is_dir() and (path/'SOURCE_MANIFEST.json').exists() and hashlib.sha256((path/'SOURCE_MANIFEST.json').read_bytes()).hexdigest()==protocol['source_bundle_sha256'])
status=dict(stage='A2_GBM',work_package='A',status='completed_pending_parent_validation' if summary_ready else 'running',
    updated_at=datetime.now(timezone.utc).isoformat(),jobs=registry,
    summary='GBM cellwise seed-mechanism replacement; original marker and DL inputs preserved. No PTC computation.',
    completed=len(units),remaining=242-len(units),unit='sample×HVG evaluation units',
    details='Counts are completed evaluation checkpoints with manifest hash checks; full scientific and work-package A acceptance remain pending.',
    completed_evaluation_units=len(units),expected_evaluation_units=242,
    expected_candidate_conditions=1210,completed_units=units,
    invalid_completion_markers=invalid_markers,
    original_regression_manifest_present=manifest_complete(OUT/'regression'),
    aggregate_manifest_present=summary_ready,whole_work_package_A_complete=False,
    evidence=[str(path.relative_to(ROOT)) for path in [OUT/'protocol.json',active_source/'PROTOCOL.md',active_source/'SOURCE_MANIFEST.json',OUT/'jobs.json',OUT/'pilot_gate.json',OUT/'initial_large_validation.json',OUT/'NL090_large_validation.json',OUT/'SN040_large_validation.json',OUT/'RESOURCE_UPDATE.json',OUT/'PROVENANCE_CLARIFICATION.md',OUT/'summary/manifest.json',OUT/'summary/REPORT.md',OUT/'validation.json'] if path.exists()],
    limitations=['Explicit seed construction replacement, not pure deletion of clustering','All five lambda values retained','Patient-based retrospective label-heldout selection','No-op and unavailable trained refinement reported explicitly'])
path=OUT/'status.json';temp=path.with_suffix('.tmp');temp.write_text(json.dumps(status,ensure_ascii=False,indent=2)+'\n');temp.replace(path)
print(json.dumps(dict(status=status['status'],completed_evaluation_units=len(units),jobs=len(registry))))
