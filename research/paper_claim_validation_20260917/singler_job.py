"""Use isolated installed SingleR environment and then the common all-cell evaluator."""
import csv
import os
from pathlib import Path
import subprocess
import sys
from common import OUT, require_slurm
if __name__=='__main__':
    require_slurm()
    import evaluate_comparator
    fold=sys.argv[1] if len(sys.argv)>1 else os.environ['SLURM_ARRAY_TASK_ID']
    source=Path(__file__).resolve().parent
    subprocess.run(['/fs/scratch/PCON0080/yimin/mamba_envs/annobench/bin/Rscript',str(source/'singler_R.R'),str(fold)],check=True)
    tests=list(csv.DictReader((OUT/'reference_inputs'/f'fold{fold}'/'test_samples.csv').open()))
    for row in tests:evaluate_comparator.run('SingleR',row['sample'])
