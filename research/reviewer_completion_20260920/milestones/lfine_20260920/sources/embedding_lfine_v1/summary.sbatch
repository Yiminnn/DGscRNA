#!/bin/bash
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=A1_Lfine_summary
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python handoff/reviewer_completion_20260920/embedding_lfine_v1/select_lfine.py
exec /fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python handoff/reviewer_completion_20260920/embedding_lfine_v1/verify_selection.py
