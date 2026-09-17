#!/bin/bash
#SBATCH --job-name=Rref_evaluate
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/Rref_evaluate_%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/r_reference_campaign_20260917/evaluate_unit.py "$@"
