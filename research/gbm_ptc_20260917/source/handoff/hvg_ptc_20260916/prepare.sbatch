#!/bin/bash
#SBATCH --job-name=gbm_hvg_prepare
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/gbm_hvg_prepare_%A_%a.out
set -euo pipefail
BASE=/fs/scratch/PCON0080/yimin/dgscrna
PY=/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 NUMBA_NUM_THREADS=8
export G274_THREADS=8 PYTHONHASHSEED=42 PYTHONOPTIMIZE=0 OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE
cd "$BASE"
"$PY" -u handoff/hvg_ptc_20260916/prepare.py --sample-list "$1" --index "$SLURM_ARRAY_TASK_ID"
