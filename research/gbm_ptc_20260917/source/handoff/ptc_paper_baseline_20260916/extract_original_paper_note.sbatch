#!/bin/bash
#SBATCH --job-name=PTC_original_note
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/PTC_original_note_%j.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/ptc_paper_baseline_20260916/extract_original_paper_note.py
