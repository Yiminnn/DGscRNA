#!/bin/bash
#SBATCH --job-name=ptc_original_workspace
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/ptc_original_workspace_%j.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript handoff/ptc_recovery_20260916/inventory_workspace.R
