#!/bin/bash
#SBATCH --job-name=PTC_original_finalize
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/PTC_original_finalize_%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/ptc_paper_baseline_20260916/finalize_full_verification.py
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/ptc_paper_baseline_20260916/correct_notebook.py
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/ptc_paper_baseline_20260916/stage_delivery.py
rclone copy results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916/GBM_PTC_results_20260916 onedrive:work_od/share/dgscrna_GSE274546_TKU3186/01_report/GBM_PTC_results_20260916 --files-from results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916/PTC_baseline_correction_upload_files.txt --transfers 4 --checkers 8 --log-file results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916/PTC_baseline_correction_final_upload.log --log-level INFO
rclone check results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916/GBM_PTC_results_20260916 onedrive:work_od/share/dgscrna_GSE274546_TKU3186/01_report/GBM_PTC_results_20260916 --files-from results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916/PTC_baseline_correction_upload_files.txt --download --one-way --checkers 4 --log-file results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916/PTC_baseline_correction_final_remote_check.log --log-level INFO
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u handoff/ptc_paper_baseline_20260916/record_final_receipt.py
