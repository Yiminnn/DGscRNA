"""Syntax validation only; no fitting is performed by this helper."""
import ast
import os
from pathlib import Path
import subprocess
from common import OUT, RSCRIPT, require_slurm, sha, write_json, utc
if __name__=='__main__':
    require_slurm();source=Path(__file__).resolve().parent
    files=[p for p in source.iterdir() if p.suffix in ['.py','.R','.sbatch']]
    for p in files:
        if p.suffix=='.py':ast.parse(p.read_text(),filename=str(p))
        if p.suffix=='.sbatch':subprocess.run(['bash','-n',str(p)],check=True)
    subprocess.run([RSCRIPT,'-e','for(p in commandArgs(trailingOnly=TRUE))parse(p)',*[str(p) for p in files if p.suffix=='.R']],check=True)
    write_json(OUT/'verification/source_syntax_checks.json',dict(status='passed',source=str(source),files={p.name:sha(p) for p in files},
        job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
