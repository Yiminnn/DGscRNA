#!/bin/bash
#SBATCH --job-name=Rref_accounting
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/Rref_accounting_%j.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/r_reference_campaign_20260917/collect_slurm_accounting.py
