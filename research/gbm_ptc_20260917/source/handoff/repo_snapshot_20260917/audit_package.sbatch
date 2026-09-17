#!/bin/bash
#SBATCH --job-name=DG_python_R_audit
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/DG_python_R_audit_%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=2
export PYTHONPATH=/fs/scratch/PCON0080/yimin/dgscrna/.worktrees/reproducibility-gbm-ptc-20260917
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/repo_snapshot_20260917/audit_package.py
cd /fs/scratch/PCON0080/yimin/dgscrna/.worktrees/reproducibility-gbm-ptc-20260917
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -m pytest tests/test_r_alignment.py -q
