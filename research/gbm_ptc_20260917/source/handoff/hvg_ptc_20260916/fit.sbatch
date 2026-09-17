#!/bin/bash
#SBATCH --job-name=gbm_hvg_fit
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/gbm_hvg_fit_%A_%a.out
set -euo pipefail
BASE=/fs/scratch/PCON0080/yimin/dgscrna
PY=/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 NUMBA_NUM_THREADS=8
export G274_THREADS=8 PYTHONHASHSEED=42 PYTHONOPTIMIZE=0 OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE
cd "$BASE"
if [ "$#" -gt 0 ]; then
  "$PY" -u handoff/hvg_ptc_20260916/fit.py "$@"
else
  "$PY" -u handoff/hvg_ptc_20260916/fit.py --task "$SLURM_ARRAY_TASK_ID"
fi
