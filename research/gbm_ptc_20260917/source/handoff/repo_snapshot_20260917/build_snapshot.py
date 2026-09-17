"""Copy the executed research source, preserving bytes and recording dependencies."""
from pathlib import Path
from datetime import datetime, timezone
import ast
import hashlib
import json
import shutil
import subprocess

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
TREE = ROOT / '.worktrees/reproducibility-gbm-ptc-20260917'
DEST = TREE / 'research/gbm_ptc_20260917'
BASE_COMMIT = '4bf17c4cb9518ab427fa92e8bf5bb4d594d30e00'
SUFFIXES = {'.py', '.R', '.r', '.rmd', '.sbatch', '.sh', '.dot'}
PRIMARY = ['hvg_ptc_20260916','gbm_cluster_figures_20260916',
           'ptc_recovery_20260916','ptc_paper_baseline_20260916']


def checksum(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main():
    paths = {p for name in PRIMARY for p in (ROOT/'handoff'/name).rglob('*')
             if p.is_file() and p.suffix in SUFFIXES and '__pycache__' not in p.parts}
    paths.update(ROOT/'handoff/ptc_recovery_20260916'/name for name in
                 ['archive_members.txt','FINAL_ANALYSIS_PLAN.md'])
    candidates = [ROOT/line for line in subprocess.check_output(
        ['rg','--files','handoff','-g','*.py'], cwd=ROOT, text=True).splitlines()]
    by_name = {}
    for path in candidates:
        by_name.setdefault(path.stem, []).append(path)
    canonical = {name:ROOT/'handoff/mc'/f'{name}.py' for name in ['bench_ext','metrics_lib','gbm_labels']}
    # Verified explicit sys.path insertion in the five PTC inventory/reconciliation scripts.
    canonical['common'] = ROOT/'handoff/hvg_ptc_20260916/common.py'
    pending = list(paths)
    processed = set()
    ambiguous = []
    imports = []
    while pending:
        path = pending.pop()
        if path in processed or path.suffix != '.py':
            continue
        processed.add(path)
        node = ast.parse(path.read_text())
        names = set()
        for part in ast.walk(node):
            if isinstance(part, ast.ImportFrom) and part.module and part.level == 0:
                names.add(part.module.split('.')[0])
            elif isinstance(part, ast.Import):
                names.update(alias.name.split('.')[0] for alias in part.names)
        for name in sorted(names):
            local = path.parent/(name+'.py')
            options = by_name.get(name, [])
            if local.exists():
                found = local
            elif name in canonical:
                found = canonical[name]
            elif len(options) == 1:
                found = options[0]
            elif len(options) > 1:
                ambiguous.append(dict(source=str(path.relative_to(ROOT)), module=name,
                                      candidates=[str(p.relative_to(ROOT)) for p in options]))
                continue
            else:
                continue
            imports.append(dict(source=str(path.relative_to(ROOT)), module=name,
                                dependency=str(found.relative_to(ROOT))))
            if found not in paths:
                paths.add(found)
                pending.append(found)
    # Include the review/probe tooling without pulling unrelated research tracks.
    paths.update(p for p in (ROOT/'handoff/repo_snapshot_20260917').iterdir()
                 if p.is_file() and p.suffix in SUFFIXES)
    DEST.mkdir(parents=True, exist_ok=True)
    records = []
    for source in sorted(paths):
        relative = source.relative_to(ROOT)
        target = DEST/'source'/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        digest = checksum(source)
        assert checksum(target) == digest
        records.append(dict(source=str(relative), snapshot=str(target.relative_to(DEST)),
                            bytes=source.stat().st_size, sha256=digest))
    legacy = ROOT/'results/hvg_ptc_20260916_v1/ptc_recovery/archive/tcr/ptc_val/scripts/DGscRNA-Share/R'
    for name in ['source.R', 'source.py']:
        source = legacy/name
        target = DEST/'reference/ptc_archive'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        records.append(dict(source=str(source.relative_to(ROOT)), snapshot=str(target.relative_to(DEST)),
                            bytes=source.stat().st_size, sha256=checksum(source)))
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(), package_base_commit=BASE_COMMIT,
                  sources_byte_preserved=True, n_files=len(records), files=records,
                  local_python_dependencies=imports, unresolved_ambiguous_imports=ambiguous,
                  excluded='Expression matrices, patient metadata, cell labels, model weights, rendered notebooks, result tables, credentials and OneDrive content',
                  portability='Executed scripts retain site-specific absolute paths. This is an archival source snapshot, not a portable standalone release.')
    (DEST/'SOURCE_MANIFEST.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['files','local_python_dependencies']}, indent=2))


if __name__ == '__main__':
    main()
