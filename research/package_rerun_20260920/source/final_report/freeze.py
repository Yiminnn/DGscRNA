"""Close source/runtime/proof metadata after independent review; no report assembly."""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import os
import sys
assert os.environ.get('SLURM_JOB_ID')
CODE = Path(__file__).resolve().parent
ROOT = CODE.parents[2]
PROOFS = ROOT / 'results/hvg_ptc_20260916_v1/package_reference_rerun_20260920/final_report_validation'
def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()
target = CODE / 'REPORT_LOCK.json'
assert not target.exists()
proofs = [PROOFS/'archived_functions_v2/verification.json', PROOFS/'guard_verification_v2.json', PROOFS/'fresh_pilot_metrics.json']
for path, status in zip(proofs, ['passed_archived_function_regression', 'passed_actual_incomplete_campaign_rejection', 'passed_fresh_pilot_metric_inputs_exact']):
    assert json.loads(path.read_text())['status'] == status
regression = json.loads(proofs[0].read_text())
for path, digest in regression['source_hashes'].items():
    assert sha(path) == digest, path
assert json.loads(proofs[1].read_text())['rejected'][0]['sha256'] == sha(CODE/'run.py')
review = json.loads((CODE/'independent_review.json').read_text())
assert review['status'] == 'passed_bounded_independent_source_review'
for path, digest in review['files'].items():
    assert sha(path) == digest, path
for path, digest in review['validation_receipts'].items():
    assert sha(path) == digest, path
paths = [p for p in CODE.iterdir() if p.is_file()] + proofs
lock = dict(status='validated_report_adapter_full_campaign_pending', files={str(p):sha(p) for p in paths},
    runtime=dict(python=str(Path(sys.executable).absolute()), versions={name:importlib.metadata.version(name)
        for name in ['numpy','pandas','scipy','matplotlib','nbformat']}),
    full_campaign_report_executed=False, job=os.environ['SLURM_JOB_ID'], step=os.environ.get('SLURM_STEP_ID'))
target.write_text(json.dumps(lock,indent=2)+'\n')
print(json.dumps(dict(status=lock['status'],lock_sha256=sha(target),files=len(paths),runtime=lock['runtime'])))
