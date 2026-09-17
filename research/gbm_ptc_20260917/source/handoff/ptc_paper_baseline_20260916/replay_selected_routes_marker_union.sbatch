#!/bin/bash
#SBATCH --job-name=PTC_original_union
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --array=0-1%2
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/PTC_original_union_%A_%a.out
set -euo pipefail
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript handoff/ptc_paper_baseline_20260916/replay_selected_routes_marker_union.R
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/ptc_paper_baseline_20260916/refine_selected_routes_marker_union.py
