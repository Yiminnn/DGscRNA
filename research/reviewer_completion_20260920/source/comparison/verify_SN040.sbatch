#!/bin/bash
#SBATCH --job-name=review_GBM_GNNcheck
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/review_GBM_GNNcheck_%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript handoff/reviewer_completion_20260920/comparison/verify_lognorm_SN040_R.R
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/reviewer_completion_20260920/comparison/verify_SN040.py
