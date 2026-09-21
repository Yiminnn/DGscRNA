#!/usr/bin/env python3
"""Finalize complete GBM reruns and their local report; never launches PTC."""
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
ROOT = HERE.parents[1]


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def need(value, message):
    if not value:
        raise RuntimeError(message)


def write(path, value):
    temporary = path.with_name(path.name + '.part.' + str(os.getpid()))
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def checked_result(path, expected_status):
    data = load(path)
    need((path.parent / 'COMPLETE').read_text().strip() == sha(path), 'Changed result completion flag')
    need(data['status'] == expected_status, 'Incomplete final result')
    for name, digest in data['outputs'].items():
        need(sha(path.parent / name) == digest, 'Changed completed result: ' + name)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', type=Path, required=True)
    parser.add_argument('--max-hours', type=float, default=46)
    args = parser.parse_args()
    need(os.environ.get('SLURM_JOB_ID'), 'Run this finalizer in SLURM')
    need(not sys.flags.optimize, 'Do not run with Python -O')
    args.campaign = args.campaign.resolve()
    sys.path.insert(0, str(HERE / 'extension_manager'))
    import common
    gate, tasks, release, root = common.gate_metadata(args.campaign)
    need(gate['scope'] == 'full_2617' and len(tasks) == 2617, 'Complete GBM scope required')
    frozen_paths = [args.campaign, Path(__file__), HERE / 'aggregate_extensions.py',
        HERE / 'embedding_summary/run.py', HERE / 'embedding_summary/numerical.py',
        HERE / 'embedding_summary/canonical_inference.py',
        HERE / 'embedding_summary/SOURCE_PROTOCOL_MANIFEST.json',
        HERE / 'embedding_summary/SUMMARY_LOCK.json',
        HERE / 'final_report/run.py', HERE / 'final_report/REPORT_LOCK.json',
        HERE / 'extension_manager/CONTROLLER_V2_LOCK.json',
        HERE / 'extension_manager/manage_v2.py', HERE / 'extension_manager/throttle_v2.py']
    frozen = {str(path): sha(path) for path in frozen_paths}
    control = root / 'control/gbm_completion'
    control.mkdir(parents=True, exist_ok=True)
    lock = control / 'LOCK'
    lock.mkdir()
    write(lock / 'owner.json', dict(job=os.environ['SLURM_JOB_ID'], pid=os.getpid(), host=os.uname().nodename))
    started = time.monotonic()
    def status(state, **fields):
        value = dict(status=state, at=datetime.now(timezone.utc).isoformat(),
            campaign_gate_sha256=frozen[str(args.campaign)], source_hashes=frozen,
            automatic_PTC=False, **fields)
        write(control / 'status.json', value)
        print(json.dumps({k:v for k,v in value.items() if k != 'source_hashes'}), flush=True)
    def verify_sources():
        for path, digest in frozen.items():
            need(sha(path) == digest, 'Finalizer source changed: ' + path)
    def execute(name, argv):
        verify_sources()
        status(name)
        with (control / (name + '.log')).open('a') as stream:
            subprocess.run(argv, cwd='/tmp', check=True, stdout=stream, stderr=subprocess.STDOUT)
    try:
        while True:
            verify_sources()
            complete_paths = [root / 'GBM_CORE_COMPLETE', root / 'GBM_EXTENSIONS_COMPLETE',
                root / 'report/notebook_receipt.json']
            if all(path.exists() for path in complete_paths):
                core, extension, report = map(load, complete_paths)
                need(core['accepted_tasks'] == 726 and extension['accepted_tasks'] == 2617,
                    'Incomplete scientific campaign')
                need(core['gate_sha256'] == gate['release_gate']['sha256']
                    and extension['campaign_gate_sha256'] == frozen[str(args.campaign)],
                    'Different scientific campaign accepted')
                need(report['status'] == 'updated' and report['core_gate_sha256'] == gate['release_gate']['sha256'],
                    'Core local notebook update not complete')
                break
            states = {}
            for name in ['core_manager', 'extension_manager', 'report_watcher']:
                path = root / 'control' / name / 'status.json'
                states[name] = load(path).get('status') if path.exists() else 'not_started'
            status('waiting_complete_science_and_core_report', stages=states)
            if time.monotonic() - started >= args.max_hours * 3600:
                status('watch_window_expired', scientific_completion_declared=False)
                return
            time.sleep(60)
        aggregate = root / 'evaluation/extensions'
        if not (aggregate / 'COMPLETE').exists():
            need(not aggregate.exists(), 'Preserve incomplete aggregation attempt for review')
            execute('aggregate_extensions', [release['runtime']['python'], '-s',
                str(HERE / 'aggregate_extensions.py'), '--campaign', str(args.campaign), '--out', str(aggregate)])
        aggregate_receipt = checked_result(aggregate / 'manifest.json', 'completed')
        need(aggregate_receipt['campaign_gate_sha256'] == frozen[str(args.campaign)]
            and aggregate_receipt['tasks'] == 2617, 'Wrong extension aggregation')
        selected = root / 'evaluation/embedding_selected'
        if not (selected / 'COMPLETE').exists():
            need(not selected.exists(), 'Preserve incomplete selection attempt for review')
            execute('select_A1_and_update_notebook', [gate['A1_runtime']['python'], '-s',
                str(HERE / 'embedding_summary/run.py'), '--campaign', str(args.campaign),
                '--aggregate', str(aggregate), '--update-notebook'])
        selected_receipt = checked_result(selected / 'manifest.json', 'passed_independent_selection')
        selected_notebook = load(selected / 'notebook_receipt.json')
        need(selected_notebook['status'] == 'updated' and selected_receipt['n_representations'] == 1694
            and selected_receipt['n_selection_rows'] == 420 and selected_receipt['n_paired_contrasts'] == 84,
            'Incomplete A1 selection/report')
        notebook = ROOT / 'notebooks/dgscrna_results.ipynb'
        final_report = root / 'evaluation/final_report'
        if not (final_report / 'COMPLETE').exists():
            need(not final_report.exists(), 'Preserve incomplete final report attempt for review')
            need(sha(notebook) == selected_notebook['notebook_candidate_sha256'],
                'Notebook changed after A1 update; review the newer edit')
            execute('revalidate_existing_reviewer_sections', [gate['A1_runtime']['python'], '-s',
                str(HERE / 'final_report/run.py'), '--campaign', str(args.campaign),
                '--aggregate', str(aggregate), '--selected', str(selected), '--update-notebook'])
        final_receipt = checked_result(final_report / 'manifest.json', 'passed_existing_sections_revalidated')
        final_notebook = load(final_report / 'notebook_receipt.json')
        need(final_notebook['status'] == 'updated'
            and sha(notebook) == final_notebook['notebook_candidate_sha256'],
            'Notebook differs from the validated final reviewer-section update')
        notebook_metadata = load(notebook)['metadata']
        reported_core = notebook_metadata.get('packaged_core_rerun', {})
        reported_A1 = notebook_metadata.get('packaged_A1_selection', {})
        reported_extensions = notebook_metadata.get('packaged_extension_report', {})
        need(reported_core.get('gate_sha256') == gate['release_gate']['sha256']
            and reported_core.get('evaluation_manifest_sha256') == sha(root / 'evaluation/manifest.json'),
            'The final notebook lost its validated core update')
        need(reported_A1.get('summary_sha256') == sha(selected / 'manifest.json'),
            'The final notebook does not contain the validated A1 update')
        need(reported_extensions.get('summary_sha256') == sha(final_report / 'manifest.json')
            and final_receipt['campaign_gate_sha256'] == frozen[str(args.campaign)],
            'The final notebook does not contain the validated extension report')
        verify_sources()
        # Bind the immutable validated candidate, so later authorized PTC
        # notebook updates do not invalidate the completed GBM evidence.
        notebook_snapshot = final_report / 'dgscrna_results_candidate.ipynb'
        need(sha(notebook_snapshot) == sha(notebook), 'GBM notebook snapshot differs')
        paths = complete_paths + [aggregate / 'manifest.json', selected / 'manifest.json',
            selected / 'notebook_receipt.json', final_report / 'manifest.json',
            final_report / 'notebook_receipt.json', notebook_snapshot]
        result = dict(status='GBM_packaged_reruns_and_Lfine_report_completed',
            completed_at=datetime.now(timezone.utc).isoformat(), job=os.environ['SLURM_JOB_ID'],
            core_tasks=726, extension_tasks=2617, core_terminal_conditions=139392,
            core_Lfine_threshold_rows=278784, extension_Lfine_threshold_rows=aggregate_receipt['threshold_rows'],
            campaign_gate_sha256=frozen[str(args.campaign)], release_gate_sha256=gate['release_gate']['sha256'],
            canonical_notebook=str(notebook), notebook_at_completion_sha256=sha(notebook),
            source_hashes=frozen, completion_hashes={str(path):sha(path) for path in paths},
            comparator_policy='Non-DG comparator predictions were not refitted by this campaign; their existing report sections are retained',
            PTC_started=False, HTML_created=False, OneDrive_accessed=False)
        target = root / 'GBM_PACKAGE_COMPLETE'
        need(not target.exists(), 'Existing completion receipt requires review')
        write(target, result)
        status('completed', receipt=str(target))
    except BaseException as error:
        status('failed_preserved', error=repr(error), automatic_retry=False)
        raise
    finally:
        (lock / 'owner.json').unlink()
        lock.rmdir()


if __name__ == '__main__':
    main()
