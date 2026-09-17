#!/bin/bash
#SBATCH --job-name=Rref_dispatch
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/Rref_dispatch_%j.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/r_reference_campaign_20260917/orchestrate.py --loop
