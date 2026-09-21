"""Reject the actual unfinished full campaign without loading scientific tables."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
assert os.environ.get('SLURM_JOB_ID')
CODE = Path(__file__).resolve().parent
ROOT = CODE.parents[2]
CAMP = ROOT / 'results/hvg_ptc_20260916_v1/package_reference_rerun_20260920'
GATE = CODE.parent / 'EXTENSION_GATE_2617.json'
assert GATE.is_file() and not (CAMP / 'GBM_EXTENSIONS_COMPLETE').exists()
def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()
notebook = ROOT / 'notebooks/dgscrna_results.ipynb'
before = sha(notebook)
commands = [
    [str(CODE / 'run.py'), '--campaign', str(GATE), '--aggregate', str(CAMP / 'evaluation/extensions'),
        '--selected', str(CAMP / 'evaluation/embedding_selected'), '--update-notebook'],
    [str(CODE.parent / 'aggregate_extensions.py'), '--campaign', str(GATE),
        '--out', str(CAMP / 'final_report_validation/forbidden_aggregate')],
]
records = []
for command in commands:
    result = subprocess.run([sys.executable, '-s'] + command, text=True, capture_output=True)
    assert result.returncode != 0 and 'Full GBM extension completion is pending' in result.stderr, result.stderr
    records.append(dict(script=command[0], sha256=sha(command[0]), returncode=result.returncode,
        error=result.stderr.splitlines()[-1]))
assert sha(notebook) == before
assert not (CAMP / 'evaluation/final_report').exists()
assert not (CAMP / 'final_report_validation/forbidden_aggregate').exists()
import run
import argparse
assert 'pandas' not in sys.modules and 'numpy' not in sys.modules
try:
    run.prerequisites(argparse.Namespace(campaign=GATE, aggregate=CAMP/'evaluation/extensions', selected=CAMP/'evaluation/embedding_selected'))
except RuntimeError as error:
    assert str(error) == 'Full GBM extension completion is pending'
else:
    raise RuntimeError('Incomplete campaign accepted')
assert 'pandas' not in sys.modules and 'numpy' not in sys.modules
value = dict(status='passed_actual_incomplete_campaign_rejection', campaign_gate_sha256=sha(GATE),
    scientific_modules_not_imported=True, notebook_unchanged_sha256=before, output_directories_not_created=True,
    rejected=records, job=os.environ['SLURM_JOB_ID'], step=os.environ.get('SLURM_STEP_ID'))
out = CAMP / 'final_report_validation'; out.mkdir(exist_ok=True)
target = out / (sys.argv[1] if len(sys.argv) == 2 else 'guard_verification.json')
assert not target.exists()
target.write_text(json.dumps(value, indent=2) + '\n')
print(json.dumps(value, indent=2))
