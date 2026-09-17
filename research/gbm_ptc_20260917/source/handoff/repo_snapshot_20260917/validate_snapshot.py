"""Verify archival source identity and syntax without importing research modules."""
from pathlib import Path
from datetime import datetime, timezone
import ast
import hashlib
import json
import os
import re
import subprocess

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
TREE = ROOT/'.worktrees/reproducibility-gbm-ptc-20260917'
SNAP = TREE/'research/gbm_ptc_20260917'
OUT = ROOT/'results/hvg_ptc_20260916_v1/python_R_audit_20260917'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    assert os.environ.get('SLURM_JOB_ID')
    manifest = json.loads((SNAP/'SOURCE_MANIFEST.json').read_text())
    assert not manifest['unresolved_ambiguous_imports']
    n_py = 0
    r_sources = []
    shell_sources = []
    for record in manifest['files']:
        source = ROOT/record['source']
        target = SNAP/record['snapshot']
        assert sha(source) == sha(target) == record['sha256'], record['snapshot']
        if target.suffix == '.py':
            ast.parse(target.read_text())
            n_py += 1
        elif target.suffix in ['.R', '.r']:
            r_sources.append(str(target))
        elif target.suffix in ['.sbatch', '.sh']:
            shell_sources.append(str(target))
    assert not subprocess.check_output(['git','diff','--name-only','HEAD'], cwd=TREE, text=True).strip(), 'An existing tracked file changed'
    for script in shell_sources:
        subprocess.run(['bash','-n',script], check=True)
    r_file = OUT/'parse_archived_sources.R'
    r_file.write_text('args <- commandArgs(trailingOnly=TRUE)\nfor (path in args) { parse(file=path) }\ncat("Parsed",length(args),"R sources\\n")\n')
    subprocess.run(['/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript',str(r_file),*r_sources], check=True)
    blocked = []
    patterns = [r'gh[pousr]_[A-Za-z0-9]{30,}', r'github_pat_[A-Za-z0-9_]{30,}',
                r'sk-[A-Za-z0-9]{30,}', r'-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----']
    for path in SNAP.rglob('*'):
        if path.is_file():
            assert path.stat().st_size < 5_000_000, f'Unexpected large artifact: {path}'
            data = path.read_text()
            if any(re.search(pattern,data) for pattern in patterns):
                blocked.append(str(path.relative_to(SNAP)))
    assert not blocked, f'Credential-like content found in {blocked}; contents withheld'
    report = dict(job=os.environ['SLURM_JOB_ID'], verified_utc=datetime.now(timezone.utc).isoformat(),
                  source_files_verified=manifest['n_files'], python_sources_parsed=n_py,
                  R_sources_parsed=len(r_sources), shell_sources_parsed=len(shell_sources),
                  all_source_hashes_match=True, existing_package_files_modified=False,
                  credential_literal_scan_passed=True,
                  scientific_probe_job='7340233', scientific_probe_job_state='COMPLETED',
                  scientific_probe_job_exit='0:0', existing_R_alignment_tests='6 passed',
                  scope='Source identity/syntax validation and previously completed synthetic data-flow probes; not a new cohort benchmark')
    (OUT/'snapshot_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__ == '__main__':
    main()
