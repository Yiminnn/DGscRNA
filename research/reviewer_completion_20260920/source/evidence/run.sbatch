#!/bin/bash
#SBATCH --job-name=review_GBM_DF
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/review_GBM_DF_%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /fs/scratch/PCON0080/yimin/dgscrna
exec /fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/reviewer_completion_20260920/evidence/gbm_evidence.py
