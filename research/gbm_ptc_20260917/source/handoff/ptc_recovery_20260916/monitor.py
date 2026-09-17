#!/usr/bin/env python3
"""Persist SLURM job and memory observations for the PTC recovery ledger."""
import json
import os
from pathlib import Path
import subprocess
import time
from datetime import datetime, timezone

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
BASE = ROOT / 'results/hvg_ptc_20260916_v1/ptc_recovery'


def run():
    assert os.environ.get('SLURM_JOB_ID')
    previous = None
    while not (BASE / 'MONITOR_STOP').exists():
        ledger = json.loads((BASE / 'job_ledger.json').read_text())
        ids = sorted({str(job) for k, v in ledger['jobs'].items()
                      if k != 'monitor' for job in [v.get('job',''), *v.get('prior_jobs',[])]
                      if str(job).isdigit()})
        if ids:
            proc = subprocess.run(['sacct', '-j', ','.join(ids), '--starttime',
                '2026-09-16T00:00:00', '-n', '-P', '-o',
                'JobID,JobName,State,ExitCode,Elapsed,ReqMem,MaxRSS,MaxVMSize,NodeList'],
                text=True, capture_output=True)
            rows = proc.stdout.strip().splitlines()
            timestamp = datetime.now(timezone.utc).isoformat()
            item = dict(timestamp=timestamp, sacct_returncode=proc.returncode,
                        fields=['JobID','JobName','State','ExitCode','Elapsed','ReqMem',
                                'MaxRSS','MaxVMSize','NodeList'],
                        rows=[r.split('|') for r in rows], stderr=proc.stderr)
            running = [r.split('|')[0] for r in rows
                       if r.split('|')[0].endswith('.batch') and r.split('|')[2] == 'RUNNING']
            if running:
                stats = subprocess.run(['sstat', '-j', ','.join(running), '-n', '-P',
                    '--format=JobID,AveRSS,MaxRSS,AveCPU,MaxDiskRead,MaxDiskWrite'],
                    text=True, capture_output=True)
                item['sstat'] = dict(fields=['JobID','AveRSS','MaxRSS','AveCPU',
                    'MaxDiskRead','MaxDiskWrite'], returncode=stats.returncode,
                    rows=[r.split('|') for r in stats.stdout.strip().splitlines()],
                    stderr=stats.stderr)
            with (BASE / 'resource_observations.jsonl').open('a') as f:
                f.write(json.dumps(item)+'\n')
            tmp = BASE / 'resource_status.json.part'
            tmp.write_text(json.dumps(item, indent=2)+'\n')
            tmp.replace(BASE / 'resource_status.json')
            state = [(r.split('|')[0], r.split('|')[2]) for r in rows
                     if '.' not in r.split('|')[0]]
            if state != previous:
                print(timestamp, state, flush=True)
                previous = state
        time.sleep(30)
    print('PTC monitor stop flag received.', flush=True)


if __name__ == '__main__':
    run()
