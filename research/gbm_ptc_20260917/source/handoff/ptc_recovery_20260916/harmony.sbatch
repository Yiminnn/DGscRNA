#!/bin/bash
#SBATCH --job-name=ptc_harmony
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --array=0-1%2
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/ptc_harmony_%A_%a.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python "$PTC_SOURCE/harmony.py"
