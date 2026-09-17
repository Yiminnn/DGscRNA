#!/bin/bash
#SBATCH --job-name=PTC_literal_DL
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --array=0-1%2
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/PTC_literal_DL_%A_%a.out
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH=/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/ptc_paper_baseline/legacy_dependencies
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/ptc_paper_baseline_20260916/verify_literal_DL.py
