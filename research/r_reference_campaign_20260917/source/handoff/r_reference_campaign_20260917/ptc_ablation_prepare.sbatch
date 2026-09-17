#!/bin/bash
#SBATCH --job-name=Rref_PTC_ablation
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=18:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/Rref_PTC_ablation_%A_%a.out
set -euo pipefail
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript handoff/r_reference_campaign_20260917/ptc_ablation_prepare.R
