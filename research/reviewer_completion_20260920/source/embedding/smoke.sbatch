#!/bin/bash
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --job-name=reviewA_geometry_gate
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /fs/scratch/PCON0080/yimin/dgscrna
exec /fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python "$1" "$2" "$3"
