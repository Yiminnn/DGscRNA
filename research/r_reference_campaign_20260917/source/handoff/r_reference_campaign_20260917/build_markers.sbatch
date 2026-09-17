#!/bin/bash
#SBATCH --job-name=Rref_markers
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/Rref_markers_%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/r_reference_campaign_20260917/build_markers.py
