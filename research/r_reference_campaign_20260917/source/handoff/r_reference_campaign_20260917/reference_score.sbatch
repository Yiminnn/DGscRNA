#!/bin/bash
#SBATCH --job-name=Rref_score
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/Rref_score_%A_%a.out
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript handoff/r_reference_campaign_20260917/reference_score.R "$@"
