#!/bin/bash
#SBATCH --job-name=ptc_R_group_cluster
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --array=0-5%6
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/ptc_R_group_cluster_%A_%a.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export PTC_CLUSTER_SOURCE="$PTC_SOURCE/cluster_groups.R"
/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript "$PTC_CLUSTER_SOURCE"
