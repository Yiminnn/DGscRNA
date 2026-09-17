#!/bin/bash
#SBATCH --job-name=ptc_independent_evaluation
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/ptc_evaluate_%A_%a.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python "$PTC_EVALUATION_SOURCE/evaluate_conditions.py"
