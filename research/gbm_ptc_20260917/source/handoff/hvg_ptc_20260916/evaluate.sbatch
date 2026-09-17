#!/bin/bash
#SBATCH --job-name=gbm_hvg_eval
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/gbm_hvg_eval_%A_%a.out
set -euo pipefail
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 NUMBA_NUM_THREADS=4 PYTHONOPTIMIZE=0
cd /fs/scratch/PCON0080/yimin/dgscrna
PY=/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python
if [ "$#" -gt 0 ]; then
  "$PY" -u handoff/hvg_ptc_20260916/evaluate.py "$@"
else
  "$PY" -u handoff/hvg_ptc_20260916/evaluate.py --index "$SLURM_ARRAY_TASK_ID"
fi
