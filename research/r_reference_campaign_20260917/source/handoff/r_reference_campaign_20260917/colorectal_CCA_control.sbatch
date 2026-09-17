#!/bin/bash
#SBATCH --job-name=Rref_CRC_CCA_control
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/Rref_CRC_CCA_control_%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DGSCRNA_ANCHOR_WORKERS=4
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript handoff/r_reference_campaign_20260917/reference_prepare.R colorectal_CCA28_ge31 colorectal
