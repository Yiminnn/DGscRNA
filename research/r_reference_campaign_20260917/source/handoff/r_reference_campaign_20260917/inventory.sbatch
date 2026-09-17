#!/bin/bash
#SBATCH --job-name=Rref_inventory
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/Rref_inventory_%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /fs/scratch/PCON0080/yimin/dgscrna
out=results/hvg_ptc_20260916_v1/r_reference_campaign_20260917
mkdir -p "$out"
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python /users/PCON0080/yimin/.agents/skills/get-available-resources/scripts/detect_resources.py -o "$out/allocation_resources.json"
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/r_reference_campaign_20260917/inventory.py
