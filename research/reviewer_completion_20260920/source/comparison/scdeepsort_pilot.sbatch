#!/bin/bash
#SBATCH --job-name=review_GBM_GNNlog
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/review_GBM_GNNlog_%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1 DGLBACKEND=pytorch
cd /fs/scratch/PCON0080/yimin/dgscrna
exec /fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/reviewer_completion_20260920/comparison/scdeepsort_corrected.py TKU4163 NL022
