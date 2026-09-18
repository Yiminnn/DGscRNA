#!/bin/bash
#SBATCH --job-name=claim_GBM
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/claim_%x_%A_%a.out
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export DGSCRNA_DEG_WORKERS=4
export PATH=/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin:$PATH
cd /fs/scratch/PCON0080/yimin/dgscrna
case "$1" in
  *.R) /fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript "$@" ;;
  *.py) /fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u "$@" ;;
  *) echo 'Expected an immutable R or Python source path' >&2; exit 2 ;;
esac
