#!/bin/bash
#SBATCH --job-name=ptc_recovery_monitor
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=512M
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/ptc_recovery_monitor_%j.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/ptc_recovery_20260916/monitor.py
