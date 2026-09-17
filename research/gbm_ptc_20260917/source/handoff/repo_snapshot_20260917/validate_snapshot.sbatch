#!/bin/bash
#SBATCH --job-name=DG_snapshot_verify
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/DG_snapshot_verify_%j.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
python handoff/repo_snapshot_20260917/validate_snapshot.py
