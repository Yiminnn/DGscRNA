"""Allow a bounded backfill window for the three pending 120k resource runs.

The original48h request remains the upper limit. Slurm may select a window of
at least36h, which exceeds1.2*completed50k_elapsed*(120000/50000)^2. This is
scheduler headroom, not an empirical runtime claim. No fitting parameters change.
"""
import json
import os
import subprocess

from common import OUT, require_slurm, checked, sha, write_json, utc
from adjust_pending_120k_memory import JOBS, inspect


def run():
    require_slurm()
    assert subprocess.check_output(['id', '-un'], text=True).strip() == 'yimin'
    path = OUT / 'verification/pending_120k_backfill_adjustment.json'
    assert not path.exists(), 'Inspect the previous adjustment before another invocation'
    state = json.loads((OUT / 'scalability_dispatch_state.json').read_text())
    pilot = state['jobs']['DG-scRNA/50000']
    assert pilot['status'] == 'complete' and checked(OUT / 'scalability/50000')
    parent = [a for a in pilot['accounting'] if a['JobID'] == pilot['job']]
    assert len(parent) == 1 and parent[0]['State'] == 'COMPLETED'
    elapsed = int(parent[0]['ElapsedRaw'])
    nominal = elapsed * (120000 / 50000) ** 2
    assert 36 * 3600 >= 1.2 * nominal
    memory_receipt = OUT / 'verification/pending_120k_memory_adjustment/receipt.json'
    memory = json.loads(memory_receipt.read_text())
    assert memory['status'] == 'completed' and memory['jobs'] == JOBS
    record = dict(status='in_progress', jobs=JOBS, min_hours=36, original_max_hours=48,
                  measured_pilot_job=pilot['job'], measured_pilot_elapsed_seconds=elapsed,
                  nominal_quadratic_reservation_seconds=nominal,
                  minimum_headroom_ratio=(36 * 3600) / nominal,
                  memory_adjustment_receipt_sha256=sha(memory_receipt), adjustments=[],
                  rationale='Three48h/384GiB jobs forecast21September start. Permit >=36h backfill with >20percent over the50k-based quadratic reservation. Not a measured complexity/runtime conclusion.',
                  source_sha256=sha(__file__), job=os.environ['SLURM_JOB_ID'], started_at=utc())
    for key, job in JOBS.items():
        assert state['jobs'][key]['job'] == job
        _, fields = inspect(job)
        assert fields['JobState'] in ['PENDING', 'RUNNING'], fields
        assert fields['MinMemoryNode'] == '384G', fields
        assert fields['TimeLimit'] == '2-00:00:00', fields
    write_json(path, record)
    for job in JOBS.values():
        before, fields = inspect(job)
        if fields['JobState'] == 'RUNNING':
            record['adjustments'].append(dict(job=job, action='unchanged_already_running', before=before))
            write_json(path, record)
            continue
        assert fields['JobState'] == 'PENDING', fields
        command = ['scontrol', 'update', f'JobId={job}', 'TimeMin=36:00:00']
        result = subprocess.run(command, text=True, capture_output=True)
        after, fields = inspect(job)
        record['adjustments'].append(dict(job=job, action='set_backfill_minimum', before=before,
                                         command=command, returncode=result.returncode,
                                         stdout=result.stdout, stderr=result.stderr, after=after,
                                         checked_at=utc()))
        write_json(path, record)
        assert result.returncode == 0, result.stderr
        assert fields['TimeMin'] == '1-12:00:00', fields
        assert fields['MinMemoryNode'] == '384G', fields
        print('PENDING_BACKFILL_MINIMUM_SET', job, fields['TimeMin'], fields['TimeLimit'], flush=True)
    record.update(status='completed', completed_at=utc())
    write_json(path, record)


if __name__ == '__main__':
    run()
