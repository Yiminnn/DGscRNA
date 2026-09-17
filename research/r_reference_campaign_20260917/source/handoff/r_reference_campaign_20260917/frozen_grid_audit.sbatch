#!/bin/bash
#SBATCH --job-name=Rref_frozen_grid_audit
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/Rref_frozen_grid_audit_%j.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python handoff/r_reference_campaign_20260917/frozen_grid_audit.py "$@"
