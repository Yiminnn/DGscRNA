#!/bin/bash
#SBATCH --job-name=gbm_hvg_tsne_diag
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/gbm_hvg_tsne_diag_%A_%a.out
set -euo pipefail
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 NUMBA_NUM_THREADS=8 PYTHONOPTIMIZE=0
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/hvg_ptc_20260916/diagnose_tsne.py --torch-runtime "$SLURM_ARRAY_TASK_ID"
