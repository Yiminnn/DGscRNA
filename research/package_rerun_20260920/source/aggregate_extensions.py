#!/usr/bin/env python3
"""Collect completed packaged control metrics without choosing configurations.

All candidate conditions remain visible. Descriptive averages do not substitute
for the separate frozen patient-fold selection used by the A1 comparison.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def need(value, message):
    if not value:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    need(os.environ.get('SLURM_JOB_ID'), 'Scientific aggregation requires SLURM')
    need(not sys.flags.optimize, 'Do not run with Python -O')
    sys.path.insert(0, str(HERE / 'extension_manager'))
    import common
    gate, tasks, release, root = common.gate_metadata(args.campaign.resolve())
    gate_hash = sha(args.campaign)
    completion_path = root / 'GBM_EXTENSIONS_COMPLETE'
    need(completion_path.exists(), 'Full GBM extension completion is pending')
    completion = load(completion_path)
    need(completion['campaign_gate_sha256'] == gate_hash,
        'Extension completion belongs to a different campaign')
    need(len(tasks) == 2617, 'The complete experiment roster is required')
    need(completion['accepted_tasks'] == len(tasks)
        and set(completion['acceptances']) == {str(t['index']) for t in tasks},
        'The complete worker acceptance roster is required')
    need(not args.out.exists(), 'Preserve an existing aggregation attempt')

    import pandas as pd
    spec = importlib.util.spec_from_file_location('frozen_core_evaluation', HERE / 'evaluate_core_lfine.py')
    provider = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(provider)
    semantics = provider.initialize()
    frames, audit, inputs = [], [], {}
    for task in tasks:
        directory = common.worker_receipt_root(root, task['index'])
        path = directory / 'acceptance.json'
        accepted = common.receipt(directory, 'acceptance.json', 'EXTENSION_VERIFIED_COMPLETE')
        need(accepted['campaign_gate_sha256'] == gate_hash and accepted['task_index'] == task['index']
            and accepted['family'] == task['family'] and accepted['sample'] == task['sample']
            and accepted['budget'] == task['budget'], 'Invalid worker acceptance')
        inputs[str(path)] = sha(path)
        need(completion['acceptances'][str(task['index'])]['sha256'] == sha(path),
            'Accepted worker receipt changed after campaign completion')
        artifact_hashes = accepted['result']['artifacts']
        metric_paths = [Path(p) for p in artifact_hashes if p.endswith('/evaluation/metrics.csv.gz')]
        need(len(metric_paths) == len(accepted['result']['evaluations']), 'Incomplete evaluation files')
        evaluation_root = Path(task.get('verification_output', task['output'])).resolve()
        need(evaluation_root.is_relative_to(root / 'extensions'), 'Unsafe frozen evaluation root')
        if 'verification_output' in task:
            need(task['family'].lower() == 'a1'
                and evaluation_root.is_relative_to(root / 'extensions/A1_verification'),
                'Only the A1 adapter uses a separate frozen verification root')
        unit_rows = 0
        for metric_path in metric_paths:
            need(metric_path.resolve().is_relative_to(evaluation_root),
                'Evaluation file escapes its task')
            need(sha(metric_path) == artifact_hashes[str(metric_path)], 'Lfine metrics changed')
            manifest_path = metric_path.parent / 'manifest.json'
            manifest = common.receipt(metric_path.parent, 'manifest.json', 'COMPLETE')
            need(sha(manifest_path) == artifact_hashes[str(manifest_path)], 'Evaluation receipt changed')
            need(manifest['frozen_semantics'] == semantics and manifest['status'] == 'completed'
                and manifest['n_valid'] == manifest['n_threshold_rows'], 'Incomplete or changed Lfine evaluation')
            frame = pd.read_csv(metric_path, dtype={'sample': str, 'patient': str, 'cutoff': str})
            need(len(frame) == manifest['n_threshold_rows'] and frame.terminal_valid.eq(True).all()
                and frame.status.eq('completed').all(), 'Incomplete Lfine rows')
            need(frame['sample'].eq(task['sample']).all() and frame.budget.eq(task['budget']).all(),
                'Evaluation belongs to another sample/budget')
            need(set(frame.stage) == {'terminal070', 'terminal090'}, 'Wrong annotation endpoints')
            frame['campaign_task'] = task['index']
            frame['campaign_family'] = task['family']
            frame['evaluation_manifest_sha256'] = sha(manifest_path)
            frames.append(frame)
            inputs[str(manifest_path)] = sha(manifest_path)
            inputs[str(metric_path)] = artifact_hashes[str(metric_path)]
            unit_rows += len(frame)
        need(unit_rows == accepted['result']['lfine_rows'], 'Task threshold-row count differs')
        audit.append(dict(index=task['index'], family=task['family'], sample=task['sample'],
            budget=task['budget'], terminal_conditions=accepted['result']['terminal_conditions'],
            lfine_rows=unit_rows, acceptance_sha256=sha(path)))

    data = pd.concat(frames, ignore_index=True)
    keys = ['campaign_family', 'sample', 'budget', 'configuration', 'route', 'library', 'cutoff', 'stage']
    need(not data.duplicated(keys).any(), 'Duplicate scientific conditions')
    need(len(data) == sum(item['lfine_rows'] for item in audit), 'Combined row count differs')
    need(len(data) == completion['lfine_rows']
        and sum(item['terminal_conditions'] for item in audit) == completion['terminal_conditions'],
        'Campaign completion totals differ from evaluated rows')
    grouping = ['campaign_family', 'family', 'budget', 'configuration', 'route', 'library', 'cutoff', 'stage']
    summary = provider.summarize(data, grouping)
    need(summary.n_samples_unavailable.eq(0).all(), 'Unavailable rows in completed experiment family')
    args.out.mkdir(parents=True)
    data.to_csv(args.out / 'metrics_all_extensions.csv.gz', index=False)
    summary.to_csv(args.out / 'descriptive_summary.csv.gz', index=False)
    pd.DataFrame(audit).to_csv(args.out / 'task_status.csv', index=False)
    for family, frame in data.groupby('campaign_family', sort=False):
        frame.to_csv(args.out / f'metrics_{family}.csv.gz', index=False)
    for path, digest in inputs.items():
        need(sha(path) == digest, 'Input changed during aggregation: ' + path)
    need(sha(args.campaign) == gate_hash, 'Campaign gate changed during aggregation')
    value = dict(status='completed', at=datetime.now(timezone.utc).isoformat(),
        job=os.environ['SLURM_JOB_ID'], campaign_gate_sha256=gate_hash,
        completion_sha256=sha(completion_path), tasks=len(tasks), threshold_rows=len(data),
        terminal_conditions=sum(item['terminal_conditions'] for item in audit),
        families={family: int(n) for family, n in data.groupby('campaign_family').size().items()},
        endpoint='Frozen compatible-target-set Lfine; terminal070 and terminal090',
        aggregation='Within-patient sample mean, then equal-weight patient mean; fixed 97/24 cohort split',
        selection_performed=False, optimality_claim=False, unavailable_conditions_omitted=False,
        family_denominators='Each frozen family has its own prespecified sample/configuration roster',
        script_sha256=sha(__file__), frozen_semantics=semantics, input_hashes=inputs,
        outputs={p.name: sha(p) for p in args.out.iterdir() if p.is_file()})
    (args.out / 'manifest.json').write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    (args.out / 'COMPLETE').write_text(sha(args.out / 'manifest.json') + '\n')
    print(json.dumps({key: value[key] for key in ['status', 'tasks', 'threshold_rows', 'terminal_conditions']}))


if __name__ == '__main__':
    main()
