#!/bin/bash
#SBATCH --job-name=gbm_hvg_singleton
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/gbm_hvg_singleton_%A_%a.out
set -euo pipefail
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 NUMBA_NUM_THREADS=8 PYTHONOPTIMIZE=0
cd /fs/scratch/PCON0080/yimin/dgscrna
PY=/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python
if [ "$#" -gt 0 ]; then
 "$PY" -u handoff/hvg_ptc_20260916/singleton_robustness.py "$@"
else
 "$PY" -u handoff/hvg_ptc_20260916/singleton_robustness.py --index "$SLURM_ARRAY_TASK_ID"
fi
