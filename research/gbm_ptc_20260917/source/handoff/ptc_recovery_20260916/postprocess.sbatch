#!/bin/bash
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/ptc_postprocess_%j.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 MPLBACKEND=Agg
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python "$PTC_POSTPROCESS_SCRIPT"
