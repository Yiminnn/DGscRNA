#!/bin/bash
#SBATCH --job-name=reviewB_GBM
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/reviewB_%A_%a.out
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1 DGSCRNA_DEG_WORKERS=2
export PYTHONPATH=/fs/scratch/PCON0080/yimin/dgscrna/handoff/paper_claim_validation_20260917
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/reviewer_completion_20260920/controls/run_controls.py "$@"
