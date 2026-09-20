#!/bin/bash
#SBATCH --job-name=review_C_inventory
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/review_C_inventory_%j.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
exec /fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/reviewer_completion_20260920/comparison/audit_inventory.py
