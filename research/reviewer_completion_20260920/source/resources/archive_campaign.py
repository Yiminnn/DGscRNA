"""Copy a bounded execution-source archive; never stage, commit, or push Git."""
from pathlib import Path
import hashlib
import json
import subprocess
from datetime import datetime, timezone

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
CODE = ROOT/'handoff/reviewer_completion_20260920'
CAMP = ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
DEST = ROOT/'DGscRNA/research/reviewer_completion_20260920'
EXTENSIONS = {'.py', '.R', '.sbatch', '.sh', '.md', '.mjs'}
EXCLUDED_NAMES = {'RESUME.md', 'RUN_STATE.md', 'COMPLETION.md', 'CURRENT_TASK.md'}
SNAPSHOTS = ['embedding_v5', 'embedding_v6', 'embedding_dispatch_v1',
             'embedding_summary_v1', 'embedding_workers4_profile_v1', 'embedding_dispatch_v2']


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    repo = ROOT/'DGscRNA'
    assert subprocess.check_output(['git', 'branch', '--show-current'], cwd=repo, text=True).strip() == 'align-r-reference'
    assert not DEST.exists(), 'Preserve existing archive; review explicitly before refreshing.'
    lock = repo/'requirements-lock.txt'
    lock_digest = digest(lock.read_bytes()) if lock.exists() else None
    collected = {}

    def add(source, target, category):
        source = Path(source)
        assert source.is_file() and not source.is_symlink(), source
        assert target not in collected, target
        data = source.read_bytes()
        collected[target] = (data, dict(path=target, source=str(source.relative_to(ROOT)),
            sha256=digest(data), bytes=len(data), category=category))

    def source_tree(source, target, category):
        # Only named small source directories; never recurse through result units/cache.
        for path in sorted(source.rglob('*')):
            if path.is_file() and '__pycache__' not in path.parts and path.suffix in EXTENSIONS and path.name not in EXCLUDED_NAMES:
                add(path, str(Path(target)/path.relative_to(source)), category)

    for folder in ['embedding', 'no_clustering', 'controls', 'comparison', 'resources', 'evidence', 'manuscript']:
        source_tree(CODE/folder, 'source/'+folder, 'current_execution_source')
    for name in ['update_notebook.py', 'progress_web.py']:
        add(CODE/name, 'source/'+name, 'report_source')
    for name in ['PLAN.md', 'REVIEWER_COVERAGE.md', 'R12_AUDIT.md', 'R3_AUDIT.md', 'SCTYPE_AUDIT.md']:
        add(ROOT/'handoff/reviewer_completion_plan_20260920'/name, 'approved_plan/'+name, 'approved_plan')
    for folder, names in {
        'controls': ['SOURCE_DERIVATION.json'],
        'no_clustering': ['SOURCE_MANIFEST.json'],
        'comparison': ['scdeepsort_remaining_samples.txt'],
    }.items():
        for name in names:
            add(CODE/folder/name, 'source/'+folder+'/'+name, 'source_provenance_or_frozen_tasks')
    for version in SNAPSHOTS:
        source = CAMP/'source_snapshots'/version
        source_tree(source, 'snapshots/'+version, 'actual_A1_execution_snapshot')
        manifest_name = 'SOURCE_MANIFEST.json' if (source/'SOURCE_MANIFEST.json').exists() else 'source_manifest.json'
        add(source/manifest_name, 'snapshots/'+version+'/'+manifest_name, 'original_source_manifest')
    source = CAMP/'no_clustering/source_v2'
    source_tree(source, 'snapshots/no_clustering_source_v2', 'actual_A2_frozen_source')
    add(source/'SOURCE_MANIFEST.json', 'snapshots/no_clustering_source_v2/SOURCE_MANIFEST.json', 'original_source_manifest')
    source_tree(CODE/'controls', 'snapshots/controls_current', 'B_current_source')
    add(CODE/'controls/SOURCE_DERIVATION.json', 'snapshots/controls_current/SOURCE_DERIVATION.json', 'B_original_source_derivation')
    for name in ['embedding.json', 'embedding_selection.json', 'embedding_parity_tasks.json', 'embedding_pilot_tasks.json', 'embedding_pilot_v4_tasks.json', 'embedding_geometry_smoke_tasks.json']:
        add(CAMP/'protocol'/name, 'protocol/'+name, 'frozen_A1_protocol')
    for source, target, category in [
        (CAMP/'no_clustering/protocol.json', 'protocol/no_clustering.json', 'frozen_A2_protocol'),
        (CAMP/'no_clustering/pilot_tasks.json', 'protocol/no_clustering_pilot_tasks.json', 'frozen_A2_tasks'),
        (CAMP/'no_clustering/full_tasks.json', 'protocol/no_clustering_full_tasks.json', 'frozen_A2_tasks'),
        (CAMP/'controls/tasks.json', 'protocol/controls_tasks.json', 'frozen_B_tasks'),
        (CAMP/'controls/verification/parity.json', 'provenance/B_default_parity.json', 'B_numerical_parity_proof'),
    ]:
        add(source, target, category)
    # Verify actual numerical bundles against their own historical manifests.
    bundle_checks = {}
    for version in SNAPSHOTS:
        source = CAMP/'source_snapshots'/version
        manifest_path = source/('SOURCE_MANIFEST.json' if (source/'SOURCE_MANIFEST.json').exists() else 'source_manifest.json')
        manifest = json.loads(manifest_path.read_text())
        hashes = manifest.get('files', {name:value for name,value in manifest.items() if Path(name).suffix in EXTENSIONS})
        assert hashes, version
        for name, expected in hashes.items():
            assert digest(collected['snapshots/'+version+'/'+name][0]) == expected, (version, name)
        bundle_checks[version] = dict(source_manifest_sha256=digest(manifest_path.read_bytes()), verified_files=len(hashes))
    a2 = json.loads((CAMP/'no_clustering/source_v2/SOURCE_MANIFEST.json').read_text())
    for name, expected in a2['files'].items():
        assert digest(collected['snapshots/no_clustering_source_v2/'+name][0]) == expected, name
    assert json.loads((CAMP/'no_clustering/protocol.json').read_text())['source_bundle_sha256'] == digest((CAMP/'no_clustering/source_v2/SOURCE_MANIFEST.json').read_bytes())
    bundle_checks['no_clustering_source_v2'] = dict(source_manifest_sha256=digest((CAMP/'no_clustering/source_v2/SOURCE_MANIFEST.json').read_bytes()), verified_files=len(a2['files']))
    parity = json.loads((CAMP/'controls/verification/parity.json').read_text())
    numerical = ['prepare_neighbors.R', 'score_neighbors.R', 'refine_learning.py']
    for name in numerical:
        assert digest(collected['snapshots/controls_current/'+name][0]) == parity['source_hashes'][name], name
    bundle_checks['controls_current'] = dict(numerical_files={n:parity['source_hashes'][n] for n in numerical},
        original_parity_job=parity['job'], numerical_source_unchanged_since_parity=True,
        orchestration_note='run_controls.py and the uniquely named control_summary.py received bookkeeping/verification updates after parity; their current hashes are recorded individually.')
    readme = '''# Reviewer completion execution-source archive

Execution is in progress. This archive preserves code and frozen protocols; it does not establish that the reviewer campaign, notebook delivery, manuscript revision, or all datasets are complete. Completion claims require checked result manifests and SLURM accounting at the site paths below.

`approved_plan/PLAN.md` defines the accepted work. `source/` records the current implementation at archive time. `snapshots/embedding_v5` and `embedding_v6` preserve the actual A1 numerical bundles; the dispatcher bundle is separate, and `snapshots/embedding_summary_v1` is the frozen final cohort selector. The selector bundled in v6 is historical and is not labelled current. `snapshots/no_clustering_source_v2` is the protocol-linked frozen A2 bundle. `snapshots/controls_current` records B, whose three numerical files exactly match its default-parity proof. `SOURCE_MANIFEST.json` records each copied byte hash, original path, and bundle verification.

The reference is the original R workflow and terminal DL/refinement endpoint. A1 compares seven representations and three clusterers with prespecified training-patient K selection. A2 replaces cluster-based seed construction with a prespecified cellwise mechanism. B tests SNN/UMAP neighbors and actual pseudo-label validation histories on three size-selected GBM pilots only; it is not a cohort-wide sensitivity result. Comparator repairs and resource accounting retain their own provenance. Legal no-op, no-known-label and actual-training endpoints remain distinct. PTC execution remains subject to the GBM completion gate. Historical SignacX 0-T-cell reporting is retained as a separately identified endpoint.

All scientific computation, parsing checks, evaluation and plots run through SLURM. These scripts retain site paths under `/fs/scratch/PCON0080/yimin/dgscrna`, site environments under `/fs/scratch/PCON0080/yimin/mamba_envs`, and dependencies on prior execution archives. This is an execution archive, not a portable or clean-environment-validated release. The original package and other branches are unchanged.

Live checked results are under `results/hvg_ptc_20260916_v1/reviewer_completion_20260920/`; historical reference results are under `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/` and `r_reference_campaign_20260917/`. The existing result notebook is `notebooks/dgscrna_results.ipynb`. Rendered notebooks, matrices, predictions, per-cell labels, patient-fold tables, model weights, environments, runtime state, credentials, private server files and remote delivery scripts are excluded from this Git archive. Referenced inputs remain outside Git.

Resource tables must distinguish allocated CPU-hours from measured CPU, summed job wall time from calendar span, requested configurations from unique actual fits, and batch MaxRSS from aggregate process peak. The separate 45-run CPU scaling experiment does not substitute for full search accounting. Current resource coverage remains partial while jobs run. `embedding_workers4_profile_v1` and `embedding_dispatch_v2` preserve the later operational profile after exact numeric equivalence and whole-job memory verification; earlier two-worker jobs remain historical. The worker probe is a separate resource experiment, not an extra scientific grid condition.
'''
    collected['README.md'] = (readme.encode(), dict(path='README.md', source='generated archive documentation', sha256=digest(readme.encode()), bytes=len(readme.encode()), category='archive_documentation'))
    DEST.mkdir(parents=True)
    for target, (data, record) in collected.items():
        path = DEST/target;path.parent.mkdir(parents=True, exist_ok=True);path.write_bytes(data)
        assert digest(path.read_bytes()) == record['sha256']
    manifest = dict(status='execution_in_progress_source_archive', archived_at=datetime.now(timezone.utc).isoformat(), branch='align-r-reference',
        scope='Only research/reviewer_completion_20260920; no Git staging/commit/push or package edits.',
        excluded=['matrices', 'predictions', 'per-cell labels', 'patient-fold tables', 'models', 'notebooks', 'environments', 'runtime state', 'credentials', 'private server files', 'remote delivery scripts'],
        portable_or_clean_environment_validated=False, scientific_completion_claimed=False,
        source_bundle_verification=bundle_checks, files=[r for _,r in collected.values()])
    (DEST/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest, indent=2)+'\n')
    assert (digest(lock.read_bytes()) if lock.exists() else None) == lock_digest
    print(json.dumps(dict(archive=str(DEST), n_files=len(collected), numerical_bundles_verified=bundle_checks, requirements_lock_unchanged=True)))


if __name__ == '__main__':
    main()
