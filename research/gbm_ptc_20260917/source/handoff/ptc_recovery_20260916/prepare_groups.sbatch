#!/bin/bash
#SBATCH --job-name=ptc_R_group_prepare
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --array=0-1%2
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=192G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/ptc_R_group_prepare_%A_%a.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript "$PTC_PREP_SOURCE"
