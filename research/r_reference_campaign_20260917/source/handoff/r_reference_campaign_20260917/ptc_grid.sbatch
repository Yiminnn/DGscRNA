#!/bin/bash
#SBATCH --job-name=Rref_PTC_grid
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --array=0-3%2
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/Rref_PTC_grid_%A_%a.out
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript handoff/r_reference_campaign_20260917/ptc_grid.R
