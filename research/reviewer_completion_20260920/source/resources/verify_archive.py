"""Read-only archive hash, syntax and exclusion checks; requires SLURM."""
import ast
from collections import Counter
from pathlib import Path
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone

assert os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_STEP_ID')
ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
DEST = ROOT/'DGscRNA/research/reviewer_completion_20260920'
OUT = ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/resources'
manifest = json.loads((DEST/'SOURCE_MANIFEST.json').read_text())
counts = Counter(); rfiles = []
for record in manifest['files']:
    path = DEST/record['path']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == record['sha256'], path
    assert path.suffix in {'.py','.R','.md','.sbatch','.sh','.mjs','.json','.txt'}, path
    if path.suffix == '.py':
        ast.parse(path.read_text(), filename=str(path)); counts['python'] += 1
    elif path.suffix == '.R':
        rfiles.append(str(path)); counts['R'] += 1
    elif path.suffix in {'.sh','.sbatch'}:
        subprocess.run(['bash','-n',str(path)],check=True); counts['shell'] += 1
    elif path.suffix == '.json':
        json.loads(path.read_text()); counts['json'] += 1
    elif path.suffix == '.mjs':
        subprocess.run(['node','--check',str(path)],check=True); counts['javascript'] += 1
assert all(p.is_file() and not p.is_symlink() for p in DEST.rglob('*') if not p.is_dir())
actual = {str(p.relative_to(DEST)) for p in DEST.rglob('*') if p.is_file()}
assert actual == {r['path'] for r in manifest['files']} | {'SOURCE_MANIFEST.json'}
if rfiles:
    subprocess.run(['/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript','--vanilla','-e',
        'for (p in commandArgs(trailingOnly=TRUE)) invisible(parse(file=p))', *rfiles],check=True)
v=dict(status='passed', at=datetime.now(timezone.utc).isoformat(), job=os.environ['SLURM_JOB_ID'], step=os.environ['SLURM_STEP_ID'],
       archive=str(DEST), source_manifest_sha256=hashlib.sha256((DEST/'SOURCE_MANIFEST.json').read_bytes()).hexdigest(),
       verified_files=len(manifest['files']), syntax_checked=dict(counts), no_execution_or_package_install=True,
       clean_environment_portability_proven=False, no_git_commit_or_push=True)
(OUT/'archive_verification.json').write_text(json.dumps(v,indent=2)+'\n')
print(json.dumps(v),flush=True)
