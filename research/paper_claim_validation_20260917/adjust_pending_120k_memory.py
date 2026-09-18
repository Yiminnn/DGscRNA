"""One bounded scheduler adjustment after the completed 50k resource pilot.

Only the three specified, still-pending 120k jobs may be changed. The existing
measured-pilot estimate includes 1.5x linear RSS and three extra condensed-distance
buffers (320.567 GiB); 384 GiB retains additional headroom. Counts, model settings,
four computational workers, source snapshot, job IDs and 48h limits are unchanged.
Original dispatcher memory_GB fields describe the initial 448 GiB submission;
this separate receipt and actual sacct records describe the effective allocation.
"""
import json
import os
import re
import subprocess
from pathlib import Path

from common import OUT, require_slurm, checked, sha, write_json, utc

JOBS = {
    'DG-scRNA/120000': '7369952',
    'DG-scRNA/120000/repeat1': '7369955',
    'DG-scRNA/120000/repeat2': '7369956',
}
MEMORY_GIB = 384


def inspect(job):
    raw = subprocess.check_output(['scontrol', 'show', 'job', job, '-o'], text=True)
    fields = dict(re.findall(r'([A-Za-z0-9_/]+)=([^\s]+)', raw))
    assert fields['JobId'] == job and fields['UserId'].startswith('yimin('), fields
    return raw, fields


def run():
    require_slurm()
    assert subprocess.check_output(['id', '-un'], text=True).strip() == 'yimin'
    dest = OUT / 'verification/pending_120k_memory_adjustment'
    path = dest / 'receipt.json'
    if path.exists():
        previous = json.loads(path.read_text())
        if previous.get('status') == 'completed':
            assert previous['jobs'] == JOBS
            for job in JOBS.values():
                _, fields = inspect(job)
                assert fields['MinMemoryNode'] == '384G', fields
            print('ALREADY_ADJUSTED', flush=True)
            return
    state = json.loads((OUT / 'scalability_dispatch_state.json').read_text())
    pilot = state['jobs']['DG-scRNA/50000']
    assert pilot['status'] == 'complete' and pilot['job'] == '7364718'
    assert checked(OUT / 'scalability/50000')
    terminal = OUT / 'GBM/SCALE_50000/hvg2000/UMAP2_HDBSCAN_R/terminal/L00_mean/terminal_manifest.json'
    assert checked(terminal.parent, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
    tm = json.loads(terminal.read_text())
    assert tm['training_executed'] is True and tm['identical_result_reused'] is False
    record = dict(status='in_progress', jobs=JOBS, target_memory_GiB=MEMORY_GIB,
                  measured_pilot=pilot, pilot_terminal_sha256=sha(terminal), adjustments=[],
                  note='Initial dispatcher memory_GB remains the original448GiB request. Use this receipt and completed sacct for the effective request. No scientific parameters or artifacts changed.',
                  source_sha256=sha(__file__), job=os.environ['SLURM_JOB_ID'], started_at=utc())
    # Check the complete target set before changing any job.
    for key, job in JOBS.items():
        planned = state['jobs'][key]
        assert planned['job'] == job and planned['pilot_gate'] == 50000
        assert planned['memory_GB'] == 448 and planned['requested_walltime'] == '48:00:00'
        assert 0 < planned['uncapped_estimate_GB'] < MEMORY_GIB
        # Keep at least15% above the already conservative distance-buffer estimate.
        assert MEMORY_GIB >= 1.15 * planned['uncapped_estimate_GB']
        source = Path(planned['source'])
        assert 'export DGSCRNA_DEG_WORKERS=4' in (source / 'job.sbatch').read_text()
        raw, fields = inspect(job)
        assert fields['JobState'] == 'PENDING' and fields['MinMemoryNode'] == '448G', fields
        assert fields['TimeLimit'] == '2-00:00:00', fields
        record.setdefault('before', {})[job] = raw
        record.setdefault('planned_estimate_GiB', {})[job] = planned['uncapped_estimate_GB']
    write_json(path, record)
    for job in JOBS.values():
        _, current = inspect(job)
        assert current['JobState'] == 'PENDING', current
        command = ['scontrol', 'update', f'JobId={job}',
                   f'MinMemoryNode={MEMORY_GIB * 1024}', 'NumCPUs=4',
                   'CPUsPerTask=4', 'MinCPUsNode=4']
        result = subprocess.run(command, text=True, capture_output=True)
        after, fields = inspect(job)
        record['adjustments'].append(dict(job=job, command=command, returncode=result.returncode,
                                         stdout=result.stdout, stderr=result.stderr, after=after,
                                         checked_at=utc()))
        write_json(path, record)
        assert result.returncode == 0, result.stderr
        assert fields['MinMemoryNode'] == '384G' and fields['TimeLimit'] == '2-00:00:00', fields
        assert int(fields['NumCPUs']) <= 98, fields
        print('PENDING_MEMORY_ADJUSTED', job, fields['NumCPUs'], fields['MinMemoryNode'], flush=True)
    record.update(status='completed', completed_at=utc())
    write_json(path, record)


if __name__ == '__main__':
    run()
