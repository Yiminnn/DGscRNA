#!/usr/bin/env python3
"""Freeze adapter metadata only after the full bounded pilot is accepted."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os

CODE = Path(__file__).resolve().parent
ROOT = CODE.parents[2]
PILOT = ROOT / 'results/hvg_ptc_20260916_v1/package_reference_rerun_20260920/embedding_adapter_v1/verification_TKU4163_hvg2000'


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    assert os.environ.get('SLURM_JOB_ID'), 'Acceptance closure requires SLURM'
    target = CODE / 'ADAPTER_LOCK.json'
    assert not target.exists(), 'Preserve an existing frozen adapter'
    manifest_path = PILOT / 'manifest.json'
    assert (PILOT / 'COMPLETE').read_text().strip() == sha(manifest_path)
    proof = json.loads(manifest_path.read_text())
    assert proof['status'] == 'passed_exact' and proof['n_representations'] == 7
    assert proof['n_partitions'] == proof['n_terminal_conditions'] == 91
    assert proof['n_Lfine_threshold_rows'] == 182
    paths = [p for p in CODE.rglob('*') if p.is_file() and
        p.suffix in {'.py', '.R', '.json', '.csv', '.md'} and '__pycache__' not in p.parts]
    paths += [manifest_path, PILOT / 'R_comparison.json']
    for path, digest in proof['source_hashes'].items():
        assert sha(path) == digest, path
        paths.append(Path(path))
    for space, record in proof['evaluation_proofs'].items():
        for key in ['parity', 'specification', 'manifest']:
            path = Path(record[key + '_path'])
            assert sha(path) == record[key + '_sha256']
            paths.append(path)
        evaluation = json.loads(Path(record['manifest_path']).read_text())
        assert evaluation['n_conditions'] == 13 and evaluation['n_valid'] == evaluation['n_threshold_rows'] == 26
        assert sha(record['metrics_path']) == record['metrics_sha256']
    runtime = json.loads((CODE / 'RUNTIME.json').read_text())
    assert sha(runtime['wheel']['path']) == runtime['wheel']['sha256']
    for space in proof['spaces']:
        config_path = PILOT.parent / 'pilot/TKU4163/hvg2000' / space / 'config.json'
        config = json.loads(config_path.read_text())
        assert config['runtime']['package_versions'] == {name: entry['version'] for name, entry in runtime['packages'].items()}
        assert config['runtime']['torch_threads'] == 4
        paths.append(config_path)
    value = dict(status='validated_bounded_adapter_not_full_campaign',
        files={str(path.resolve()): sha(path) for path in sorted(set(paths))},
        pilot=dict(path=str(manifest_path), sha256=sha(manifest_path), representations=7,
            terminal_conditions=91, Lfine_threshold_rows=182),
        planned_full_campaign=dict(tasks=1694, terminal_conditions=22022, threshold_rows=44044, submitted=False),
        ICA_scope='Pilot exercises parallel5000; fallback policy source is unchanged and requires per-unit observed verification',
        definition='A1 UMAP directly uses scaled HVGs; package default UMAP uses PCA30',
        created_at=datetime.now(timezone.utc).isoformat(), job=os.environ['SLURM_JOB_ID'])
    target.write_text(json.dumps(value, indent=2) + '\n')
    print(json.dumps(dict(status=value['status'], files=len(value['files']), lock_sha256=sha(target), pilot=value['pilot'])))


if __name__ == '__main__':
    main()
