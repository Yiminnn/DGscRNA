#!/bin/bash
#SBATCH --job-name=ptc_R_DL_annotation
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/ptc_annotation_%A_%a.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4
export PTC_SCORE_SOURCE="$PTC_SCORE_SNAPSHOT/score_partitions.R"
/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript "$PTC_SCORE_SOURCE" "$PTC_MODE"
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python "$PTC_DL_SNAPSHOT/refine.py" "$PTC_MODE"
