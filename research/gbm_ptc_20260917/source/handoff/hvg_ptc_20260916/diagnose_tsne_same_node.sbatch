#!/bin/bash
#SBATCH --job-name=gbm_hvg_tsne_control
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/gbm_hvg_tsne_control_%j.out
set -euo pipefail
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 NUMBA_NUM_THREADS=8 PYTHONOPTIMIZE=0 PYTHONHASHSEED=42
cd /fs/scratch/PCON0080/yimin/dgscrna
PY=/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python
for flag in 0 1; do
 "$PY" -u handoff/hvg_ptc_20260916/diagnose_tsne.py --torch-runtime "$flag" --run-label same_node
done
