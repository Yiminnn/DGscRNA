#!/usr/bin/env python3
"""Wait for complete core acceptance, then update the local notebook in SLURM."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gate', type=Path, required=True)
    parser.add_argument('--plot-python', type=Path, required=True)
    parser.add_argument('--max-hours', type=float, default=46)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Run report watcher in SLURM')
    args.gate = args.gate.resolve(strict=True)
    # Preserve a virtual-environment launcher symlink while making cwd irrelevant.
    args.plot_python = args.plot_python.absolute()
    if not args.plot_python.is_file() or not os.access(args.plot_python, os.X_OK):
        raise RuntimeError('Plotting Python is not executable')
    frozen = {str(path): sha(path) for path in [args.gate, Path(__file__),
        HERE / 'refresh_core_notebook.py',
        HERE.parents[1] / 'handoff/lfine_compact_20260920/plot_shared_umap.py']}
    gate = json.loads(args.gate.read_text())
    root = Path(gate['output_root'])
    control = root / 'control/report_watcher'
    control.mkdir(parents=True, exist_ok=True)
    lock = control / 'LOCK'
    lock.mkdir()  # An existing controller needs explicit owner inspection.
    (lock / 'owner.json').write_text(json.dumps(dict(pid=os.getpid(),
        job=os.environ['SLURM_JOB_ID'], host=os.uname().nodename)) + '\n')
    start = time.monotonic()
    def status(state, **fields):
        value = dict(status=state, at=datetime.now(timezone.utc).isoformat(),
            frozen_sources=frozen, **fields)
        temporary = control / 'status.json.part'
        temporary.write_text(json.dumps(value, indent=2) + '\n')
        temporary.replace(control / 'status.json')
        print(json.dumps({k: v for k, v in value.items() if k != 'frozen_sources'}), flush=True)
    try:
        while True:
            if any(sha(path) != digest for path, digest in frozen.items()):
                raise RuntimeError('Frozen report source changed; explicit restart required')
            if (root / 'GBM_CORE_COMPLETE').exists():
                status('assembling_complete_core_report')
                with (control / 'refresh.log').open('a') as stream:
                    subprocess.run([str(args.plot_python), '-s', str(HERE / 'refresh_core_notebook.py'),
                        '--gate', str(args.gate.resolve())], check=True, cwd='/tmp',
                        stdout=stream, stderr=subprocess.STDOUT)
                status('completed', report_receipt=str(root / 'report/notebook_receipt.json'))
                return
            manager_path = root / 'control/core_manager/status.json'
            manager = json.loads(manager_path.read_text()) if manager_path.exists() else {}
            if manager.get('status', '').startswith(('blocked', 'failed', 'error')):
                status('waiting_for_core_repair', core_status=manager)
            else:
                status('waiting_for_all_726_core_units', verified=manager.get('verified_tasks', 0))
            if time.monotonic() - start >= args.max_hours * 3600:
                status('watch_window_expired', notebook_updated=False)
                return
            time.sleep(60)
    except BaseException as error:
        status('failed', error=repr(error))
        raise
    finally:
        (lock / 'owner.json').unlink()
        lock.rmdir()


if __name__ == '__main__':
    main()
