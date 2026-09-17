#!/bin/bash
#SBATCH --job-name=gbm_notebook_extend
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/gbm_notebook_extend_%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONHASHSEED=42 PYTHONOPTIMIZE=0 MPLBACKEND=Agg
export JUPYTER_PATH=/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/jupyter
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/gbm_cluster_figures_20260916/assemble.py
