#!/bin/bash
#SBATCH --job-name=gbm_hvg_numeric
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/gbm_hvg_numeric_%j.out
set -euo pipefail
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 NUMBA_NUM_THREADS=8 PYTHONOPTIMIZE=0
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/hvg_ptc_20260916/retry_numerical.py
