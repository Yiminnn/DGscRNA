#!/bin/bash
#SBATCH --job-name=ptc_matched_R_single
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --array=0-7%8
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/ptc_matched_R_single_%A_%a.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
export PTC_MATCHED_SOURCE="$PTC_MATCHED_SNAPSHOT/cluster_matched_samples.R"
export PTC_SCORE_SOURCE="$PTC_SCORE_SNAPSHOT/score_partitions.R"
export PTC_MATCHED_SCORE="$PTC_MATCHED_SNAPSHOT/score_matched_samples.R"
/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript "$PTC_MATCHED_SOURCE"
/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript "$PTC_MATCHED_SCORE" matched
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python "$PTC_MATCHED_SNAPSHOT/refine_matched_samples.py"
