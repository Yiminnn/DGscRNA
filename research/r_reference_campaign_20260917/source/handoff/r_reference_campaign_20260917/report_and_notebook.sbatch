#!/bin/bash
#SBATCH --job-name=Rref_report_notebook
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/Rref_report_notebook_%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PATH=/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin:$PATH
cd /fs/scratch/PCON0080/yimin/dgscrna
python -u handoff/r_reference_campaign_20260917/report_and_notebook.py "$@"
