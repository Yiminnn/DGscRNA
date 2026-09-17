"""Collect only campaign-referenced SLURM jobs and preserve failed-attempt logs."""
import csv
import datetime
import hashlib
import json
import os
from pathlib import Path
import pwd
import re
import shutil
import subprocess
import sys

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT = ROOT / 'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'
DEST = OUT / 'resources'
if __name__ == '__main__' and not os.environ.get('DGSCRNA_ACCOUNTING_EXECUTION_SOURCE'):
    assert os.environ.get('SLURM_JOB_ID')
    snapshot_dir = OUT / 'execution_sources' / os.environ['SLURM_JOB_ID']
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    snapshot = snapshot_dir / 'collect_slurm_accounting.py'
    if not snapshot.exists():
        temporary = snapshot.with_name(snapshot.name + '.part')
        shutil.copy2(__file__, temporary)
        temporary.replace(snapshot)
    os.environ['DGSCRNA_ACCOUNTING_EXECUTION_SOURCE'] = str(snapshot)
    os.execv(sys.executable, [sys.executable, str(snapshot), *sys.argv[1:]])
JOB_ID = re.compile(r'^[0-9]{5,12}(?:_[0-9]+)?$')
FIELDS = ['JobID', 'JobIDRaw', 'JobName', 'User', 'State', 'ExitCode', 'Start',
          'End', 'Elapsed', 'ElapsedRaw', 'AllocCPUS', 'ReqMem', 'MaxRSS',
          'MaxRSSNode', 'MaxRSSTask', 'NNodes', 'Partition', 'StdOut', 'StdErr']
TERMINAL_STATES = {'COMPLETED', 'CANCELLED', 'FAILED', 'TIMEOUT', 'OUT_OF_MEMORY',
                   'NODE_FAIL', 'PREEMPTED', 'BOOT_FAIL', 'DEADLINE', 'REVOKED'}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, value):
    temporary = path.with_name(path.name + '.part')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def atomic_csv(path, rows, fields):
    temporary = path.with_name(path.name + '.part')
    with temporary.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def run():
    assert os.environ.get('SLURM_JOB_ID'), 'Run accounting collection through SLURM'
    assert pwd.getpwuid(os.geteuid()).pw_name == 'yimin'
    DEST.mkdir(parents=True, exist_ok=True)
    if (DEST / 'accounting_manifest.json').exists():
        previous = json.loads((DEST / 'accounting_manifest.json').read_text())
        history = DEST / 'history' / (previous['job'] + '_' + sha(DEST / 'accounting_manifest.json')[:12])
        history.mkdir(parents=True, exist_ok=True)
        for name in ['accounting_manifest.json', *previous['files']]:
            target = history / name
            if target.exists():
                assert sha(target) == sha(DEST / name)
            else:
                shutil.copy2(DEST / name, target)
    origins = {}
    warnings = []

    def add(value, source):
        value = str(value)
        if JOB_ID.fullmatch(value):
            origins.setdefault(value, set()).add(source)

    def discover(value, source, context=False):
        if isinstance(value, dict):
            for key, child in value.items():
                is_job = 'job' in key.lower() or key.lower() in {'retry_of', 'completion_retry'}
                discover(child, source + ':' + key, context or is_job)
        elif isinstance(value, list):
            for child in value:
                discover(child, source, context)
        elif context:
            add(value, source)

    def read_metadata(path):
        try:
            discover(json.loads(path.read_text()), str(path.relative_to(OUT)))
        except json.JSONDecodeError as exc:
            # In-progress manifests may be partially written; report, never infer IDs.
            warnings.append({'path': str(path), 'error': str(exc)})

    for path in sorted(OUT.glob('*events.jsonl')):
        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            try:
                discover(json.loads(line), f'{path.name}:{lineno}')
            except json.JSONDecodeError as exc:
                warnings.append({'path': str(path), 'line': lineno, 'error': str(exc)})
    primary_files = set(OUT.glob('*.json'))
    primary_files.update((OUT / 'verification').rglob('*.json'))
    primary_files.update((OUT / 'summary').glob('*.json'))
    for path in sorted(primary_files):
        read_metadata(path)
    add(os.environ['SLURM_JOB_ID'], 'current_accounting_collection')
    primary_ids = set(origins)
    # Fitting/training manifests recover manually submitted and array-task jobs.
    for path in sorted(OUT.rglob('*manifest.json')):
        if DEST not in path.parents and path not in primary_files:
            read_metadata(path)
    for path in (OUT / 'execution_sources').iterdir():
        if path.is_dir():
            add(path.name, 'execution_sources/' + path.name)

    rows_by_id = {}
    queries = []

    def collect(job_ids):
        ordered = sorted(job_ids, key=lambda value: (int(value.split('_')[0]), value))
        for start in range(0, len(ordered), 150):
            batch = ordered[start:start + 150]
            command = ['sacct', '--noheader', '--parsable2', '--array', '--units=K', '--user=yimin',
                       '--jobs=' + ','.join(batch),
                       '--format=' + ','.join(field + '%320' for field in FIELDS)]
            result = subprocess.run(command, text=True, capture_output=True, check=True)
            queries.append({'job_ids': batch, 'command': command, 'stderr': result.stderr})
            for row in csv.DictReader(result.stdout.splitlines(), fieldnames=FIELDS, delimiter='|'):
                if not row['JobID']:
                    continue
                if row.get(None):
                    raise ValueError('Unexpected sacct columns: ' + repr(row))
                assert row['User'] in {'', 'yimin'}, 'Refusing accounting for another identity'
                rows_by_id[row['JobID']] = row

    collect(primary_ids)
    # Querying an array parent already returns its task rows and allocation IDs.
    # Avoid thousands of redundant sacct calls for those training-manifest IDs.
    covered = set()
    for row in rows_by_id.values():
        covered.update([row['JobID'].split('.')[0], row['JobIDRaw'].split('.')[0]])
    collect(set(origins) - covered - primary_ids)
    rows = sorted(rows_by_id.values(), key=lambda row: row['JobID'])
    atomic_csv(DEST / 'slurm_accounting_steps.csv', rows, FIELDS)
    steps = {}
    for row in rows:
        steps.setdefault(row['JobID'].split('.')[0], []).append(row)

    def rss_kib(value):
        if not value:
            return None
        match = re.fullmatch(r'([0-9.]+)([KMGTPE]?)(?:i?B)?', value)
        if not match:
            raise ValueError('Unrecognized MaxRSS: ' + value)
        # sacct --units=K emits K, but support the documented suffixes explicitly.
        factors = {'': 1, 'K': 1, 'M': 1024, 'G': 1024 ** 2,
                   'T': 1024 ** 3, 'P': 1024 ** 4, 'E': 1024 ** 5}
        return float(match[1]) * factors[match[2]]

    ledger = []
    failed_logs = []
    missing_logs = []
    preserved = {}
    failed_dir = DEST / 'failed_logs'
    failed_dir.mkdir(exist_ok=True)
    for row in rows:
        if '.' in row['JobID']:
            continue
        job = row['JobID']
        state = row['State'].split()[0].rstrip('+')
        final = state in TERMINAL_STATES
        task_steps = steps[job]
        measurements = [(rss_kib(item['MaxRSS']), item['JobID']) for item in task_steps if item['MaxRSS']]
        peak, peak_step = max(measurements, default=(None, ''), key=lambda item: item[0])
        failed = final and (state != 'COMPLETED' or row['ExitCode'] != '0:0')
        item = {key: row[key] for key in ['JobID', 'JobIDRaw', 'JobName', 'State', 'ExitCode',
                'Start', 'End', 'Elapsed', 'ElapsedRaw', 'AllocCPUS', 'ReqMem', 'NNodes', 'Partition']}
        item.update(accounting_final=final, accounting_status='final' if final else 'pending_or_running',
                    failed_attempt=failed, MaxRSS_KiB=peak, MaxRSS_source_step=peak_step,
                    MaxRSS_available=peak is not None)
        ledger.append(item)
        if not failed:
            continue
        parent, _, task = job.partition('_')
        candidates = set()
        for field in ['StdOut', 'StdErr']:
            pattern = row[field]
            for token, replacement in {'%A': parent, '%a': task or '4294967294',
                    '%j': row['JobIDRaw'], '%u': 'yimin', '%x': row['JobName']}.items():
                pattern = pattern.replace(token, replacement)
            if pattern and '%' not in pattern:
                candidate = Path(pattern)
                if candidate.is_relative_to(ROOT / 'logs'):
                    candidates.add(candidate)
        found = [path for path in candidates if path.is_file()]
        if not found:
            missing_logs.append({'job': job, 'state': row['State'],
                                 'reason': 'No log exists at the accounting paths; cancelled-before-start jobs may have none.',
                                 'expected_paths': [str(path) for path in sorted(candidates)]})
        for source in sorted(found):
            target = failed_dir / source.name
            source_hash = sha(source)
            if target.exists() and sha(target) != source_hash:
                # Preserve any earlier snapshot rather than replacing distinct content.
                target = failed_dir / (source.stem + '_' + source_hash[:12] + source.suffix)
            if not target.exists():
                shutil.copy2(source, target)
            assert sha(target) == source_hash
            key = str(target.relative_to(OUT))
            preserved[key] = {'path': key, 'source': str(source), 'sha256': source_hash,
                              'bytes': target.stat().st_size}
            failed_logs.append({'job': job, 'state': row['State'], **preserved[key]})
    ledger_fields = ['JobID', 'JobIDRaw', 'JobName', 'State', 'ExitCode', 'Start', 'End',
                    'Elapsed', 'ElapsedRaw', 'AllocCPUS', 'ReqMem', 'NNodes', 'Partition',
                    'accounting_final', 'accounting_status', 'failed_attempt', 'MaxRSS_KiB',
                    'MaxRSS_source_step', 'MaxRSS_available']
    atomic_csv(DEST / 'slurm_job_ledger.csv', ledger, ledger_fields)
    atomic_json(DEST / 'job_origin_index.json', {job: sorted(values) for job, values in sorted(origins.items())})
    atomic_json(DEST / 'accounting_queries.json', queries)
    atomic_json(DEST / 'failed_job_logs.json', {'logs': failed_logs, 'missing_logs': missing_logs})
    represented = {row[key].split('.')[0] for row in rows for key in ['JobID', 'JobIDRaw']}
    represented.update(row['JobID'].split('_')[0] for row in rows)
    timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    report = dict(status='collected', collected_at_utc=timestamp, job=os.environ['SLURM_JOB_ID'],
        identity='yimin', source_sha256=sha(__file__), execution_source=str(Path(__file__).resolve()),
        scope='Only job IDs recorded by campaign events, manifests, retry records, source snapshots or this collection job.',
        discovered_job_references=len(origins), allocation_or_array_task_rows=len(ledger), accounting_step_rows=len(rows),
        final_rows=sum(row['accounting_final'] for row in ledger),
        pending_or_running_rows=sum(not row['accounting_final'] for row in ledger),
        failed_attempt_rows=sum(row['failed_attempt'] for row in ledger),
        missing_accounting_ids=sorted(set(origins) - represented),
        missing_or_partial_metadata=warnings, failed_log_files=len(preserved),
        failed_log_records=len(failed_logs), failed_jobs_without_log=len(missing_logs),
        accounting_limitations=['MaxRSS comes from the maximum available task step, commonly .batch; it is not summed across steps or an aggregate node peak.',
            'Blank MaxRSS remains unavailable, never zero. Pending/running rows are explicitly non-final.',
            'Completed and failed attempts are retained separately; cancelled-before-start jobs may have no resource measurement or log.',
            'This collector and the controller/finalizer can still be active at the delivery snapshot; finalizer exit is independently verified afterward.'],
        files={name: sha(DEST / name) for name in ['slurm_accounting_steps.csv', 'slurm_job_ledger.csv',
               'job_origin_index.json', 'accounting_queries.json', 'failed_job_logs.json']}, failed_logs=list(preserved.values()))
    atomic_json(DEST / 'accounting_manifest.json', report)
    print(json.dumps({key: report[key] for key in ['status', 'allocation_or_array_task_rows', 'final_rows',
          'pending_or_running_rows', 'failed_attempt_rows', 'failed_log_files', 'missing_accounting_ids']}, indent=2), flush=True)


if __name__ == '__main__':
    run()
