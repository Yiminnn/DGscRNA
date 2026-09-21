#!/usr/bin/env python3
"""Metadata-only waiter; launches exactly one bounded scientific SLURM verifier."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import subprocess
import time

CODE = Path(__file__).resolve().parent
ROOT = CODE.parents[2]
CAMP = ROOT / 'results/hvg_ptc_20260916_v1/package_reference_rerun_20260920'
ADAPTER = CAMP / 'embedding_adapter_v1'
ACTUAL = ADAPTER / 'pilot/TKU4163/hvg2000'
EXPECTED = ROOT / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/embedding/TKU4163/hvg2000'
OUT = ADAPTER / 'verification_TKU4163_hvg2000'
PYTHON = CAMP / 'runtime_embedding_v1/bin/python'
EVAL_PYTHON = CAMP / 'runtime_v3/bin/python'
RSCRIPT = ROOT / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/reference_examples/isolated_7447099/r/bin/Rscript'
SPACES = ['noDR', 'PCA2', 'FA2', 'ICA2', 'Isomap2', 'UMAP2', 'TSNE2']


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def status(value):
    value.update(updated_at=datetime.now(timezone.utc).isoformat(), full_campaign_submitted=False)
    path = ADAPTER / 'ADAPTER_STATUS.json'
    temporary = path.with_suffix('.part.json')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def main():
    started = time.monotonic()
    while True:
        complete = []
        for space in SPACES:
            directory = ACTUAL / space
            flag, manifest = directory / 'FIT_COMPLETE', directory / 'fit_manifest.json'
            if flag.exists() and manifest.exists() and flag.read_text().strip() == sha(manifest):
                complete.append(space)
        status(dict(status='pilot_fitting', complete_spaces=complete, expected_spaces=SPACES,
            valid_fit_units=len(complete), expected_fit_units=7, accepted_terminal_conditions=0,
            note='Completed fitting units still require independent terminal verification'))
        if len(complete) == 7:
            break
        log = ROOT / 'logs/pkg_embedding_TKU4163_pilot_v2.log'
        if log.exists() and 'Traceback (most recent call last):' in log.read_text():
            status(dict(status='pilot_failed_preserved', complete_spaces=complete, log=str(log)))
            raise SystemExit(2)
        if time.monotonic() - started > 3 * 3600:
            status(dict(status='watch_timeout_fit_not_accepted', complete_spaces=complete))
            raise SystemExit(3)
        time.sleep(30)
    status(dict(status='independent_verification_running', valid_fit_units=7, expected_terminal_conditions=91))
    command = ['srun', '--jobid=7204040', '--overlap', '--exact', '--nodes=1', '--ntasks=1',
        '--cpus-per-task=2', '--mem=8G', '--time=00:30:00', str(PYTHON), '-s', str(CODE / 'verify.py'),
        '--actual', str(ACTUAL), '--expected', str(EXPECTED), '--rscript', str(RSCRIPT),
        '--evaluation-python', str(EVAL_PYTHON), '--out', str(OUT)]
    with (ROOT / 'logs/pkg_embedding_TKU4163_verification.log').open('x') as log:
        result = subprocess.run(command, cwd='/tmp', stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        status(dict(status='verification_failed_preserved', returncode=result.returncode,
            output=str(OUT), log=str(ROOT / 'logs/pkg_embedding_TKU4163_verification.log')))
        raise SystemExit(result.returncode)
    manifest = json.loads((OUT / 'manifest.json').read_text())
    assert (OUT / 'COMPLETE').read_text().strip() == sha(OUT / 'manifest.json')
    assert manifest['status'] == 'passed_exact' and manifest['n_partitions'] == 91
    status(dict(status='bounded_pilot_passed_exact', output=str(OUT),
        manifest_sha256=sha(OUT / 'manifest.json'), n_representations=7,
        n_terminal_conditions=91, n_Lfine_threshold_rows=182))
    print('BOUNDED_A1_PILOT_VERIFIED', str(OUT), flush=True)


if __name__ == '__main__':
    main()
