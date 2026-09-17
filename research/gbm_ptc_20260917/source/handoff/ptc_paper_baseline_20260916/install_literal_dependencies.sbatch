#!/bin/bash
#SBATCH --job-name=PTC_legacy_deps
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/PTC_legacy_deps_%j.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -m pip install --no-deps --target results/hvg_ptc_20260916_v1/ptc_paper_baseline/legacy_dependencies loompy==3.0.8 numpy-groupies==0.11.3
