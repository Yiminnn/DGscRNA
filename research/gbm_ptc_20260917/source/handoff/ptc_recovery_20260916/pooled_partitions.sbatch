#!/bin/bash
#SBATCH --job-name=ptc_pooled_partition
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --array=0-31%32
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/ptc_pooled_partition_%A_%a.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
export PTC_SCORE_SOURCE="$PTC_SCORE_SNAPSHOT/score_partitions.R"
export PTC_PARTITION_SCORE="$PTC_PARTITION_SNAPSHOT/score_pooled_partition.R"
/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript "$PTC_PARTITION_SCORE" partition
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python "$PTC_PARTITION_SNAPSHOT/refine_pooled_partition.py"
