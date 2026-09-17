#!/bin/bash
#SBATCH --job-name=Rref_density_parity
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/Rref_density_parity_%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript handoff/r_reference_campaign_20260917/validate_density.R
