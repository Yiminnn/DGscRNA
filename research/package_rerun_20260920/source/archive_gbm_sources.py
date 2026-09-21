#!/usr/bin/env python3
"""Archive validated rerun code in the authorized repository, without Git actions."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REPO = ROOT / 'DGscRNA'
DEST = REPO / 'research/package_rerun_20260920'


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', type=Path, required=True)
    args = parser.parse_args()
    gate = json.loads(args.campaign.read_text())
    if gate.get('status') != 'reviewed_extension_campaign' or gate.get('scope') != 'full_2617':
        raise RuntimeError('Archive only the reviewed full GBM campaign')
    if DEST.exists():
        raise RuntimeError('Preserve existing source archive')
    if not (HERE / 'final_report/REPORT_LOCK.json').is_file():
        raise RuntimeError('Finish and validate the final report adapter before archiving it')
    top = ['EXECUTION_PLAN.md', 'rerun_inventory.md', 'rerun_inventory.json',
        'rerun_inventory.core_tasks.json', 'RELEASE_GATE.json', 'core_task.py',
        'core_array.sbatch', 'manage_core.py', 'verify_packaged_parity.py',
        'verify_packaged_r_artifacts.R', 'evaluate_core_lfine.py',
        'evaluate_core_lfine_sources.json', 'evaluate_extension_lfine.py',
        'geometry_control.py', 'GEOMETRY_CONTROL.md', 'mlp_controls_adapter.py',
        'refine_controls_packaged.py', 'mlp_controls_tasks.json', 'mlp_controls_sources.json',
        'mlp_controls_array.sbatch', 'MLP_CONTROLS.md', 'representation_adapter.py',
        'representation_configurations.py', 'representation_sources.json',
        'representation_tasks.json', 'representation_derivation.json',
        'prepare_representation_R.R', 'score_representation_R.R',
        'verify_representation_inputs.R', 'verify_representation_roundtrip.R',
        'LFINE_EVALUATION.md', 'refresh_core_notebook.py', 'watch_core_report.py',
        'aggregate_extensions.py', 'watch_gbm_completion.py', 'archive_gbm_sources.py',
        'review_extension_candidate.py', 'gbm_completion_guard_review.json']
    directories = ['no_cluster_adapter', 'reviewer_b_adapter', 'embedding_adapter',
        'embedding_summary', 'extension_manager', 'final_report']
    paths = [HERE / name for name in top if (HERE / name).is_file()]
    for name in directories:
        paths += [p for p in (HERE / name).rglob('*') if p.is_file()
            and p.suffix in {'.py', '.R', '.json', '.csv', '.md', '.sbatch'}
            and '__pycache__' not in p.parts and 'original' not in p.parts
            and not p.name.startswith('CAMPAIGN_CANDIDATE')]
    # Include exact frozen numerical originals and their source manifests.
    for name in ['embedding_adapter/original', 'embedding_summary/original']:
        paths += [p for p in (HERE / name).rglob('*') if p.is_file()
            and p.suffix in {'.py', '.R', '.json'} and '__pycache__' not in p.parts]
    for path, digest in gate['source_files'].items():
        source = Path(path).resolve()
        if source.is_relative_to(HERE):
            if sha(source) != digest:
                raise RuntimeError('Frozen source changed: ' + str(source))
            paths.append(source)
    paths.append(args.campaign.resolve())
    paths.append(Path(gate['task_manifest']['path']).resolve())
    paths.append(Path(gate['root_review']['candidate']['path']).resolve())
    paths = sorted(set(paths))
    inventory = {}
    for source in paths:
        relative = source.relative_to(HERE)
        destination = DEST / 'source' / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        digest = sha(source)
        if sha(destination) != digest:
            raise RuntimeError('Copy verification failed: ' + str(source))
        inventory[str(relative)] = dict(sha256=digest, original=str(source))
    manifest = dict(status='validated_sources_archived_full_rerun_still_in_progress',
        archived_at=datetime.now(timezone.utc).isoformat(), release='v2.0.0rc1',
        release_commit='42e332f617ddb7196ba89529734ceebcdd068b9e',
        campaign_gate_sha256=sha(args.campaign), core_units=726, extension_units=2617,
        raw_matrices_copied=False, trained_models_copied=False, notebook_copied=False,
        files=inventory)
    (DEST / 'SOURCE_MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
    (DEST / 'README.md').write_text('''# Installed-package GBM rerun sources

The published reference package is `dgscrna==2.0.0rc1`:
https://github.com/Yiminnn/DGscRNA/releases/tag/v2.0.0rc1

This archive records the SLURM orchestration, independent verification, frozen
Lfine evaluation and research ablation adapters used to rerun the existing GBM
experiments through that installed release. It preserves the exact execution
sources, including recorded local paths. For a new user's analysis, use the
portable package API and `docs/reference_workflow.md` in the repository.

Scope: 726 core sample/feature units (139,392 terminal configurations), plus 2,617
geometry, MLP, representation, neighborhood, learning, cellwise-seed and seven-space
comparison units. All new terminal results use the original Python DL/refinement
after the original R scientific stages. Marker-only calls and partition metrics
are not substituted for final annotation results.

The source archive and passing pilots do **not** establish full-cohort completion
or optimality. Execution receipts under the local result directory distinguish
pending, failed, valid no-training and trained conditions. Full GBM completion
requires separate core, extension, evaluation and notebook receipts. PTC and
public reviewer reruns follow their own baseline/input contracts and are not
certified by this GBM archive.

The active extension scheduler is `manage_v2.py`. Its separately reviewed
`CONTROLLER_V2_LOCK.json` preserves the original scientific gate and workers.
It handles SLURM's nonzero reply when a throttle update succeeds but completed
array members also produce diagnostics, verifies the actual cap, and records
the preserved state migration. The original scheduler source is retained.

The two UMAP experiments differ intentionally: the original anchor uses PCA30
before UMAP; the seven-space comparison uses scaled HVGs directly. Dataset-specific
preprocessing and the original 2,000-HVG default are preserved. All-gene and other
feature budgets are comparison conditions, not a uniform-optimality claim.

`SOURCE_MANIFEST.json` records byte hashes. Expression matrices, trained weights,
full annotations and notebook binaries are not included in this source snapshot.
''')
    print(json.dumps(dict(destination=str(DEST), files=len(inventory),
        manifest_sha256=sha(DEST / 'SOURCE_MANIFEST.json'))))


if __name__ == '__main__':
    main()
