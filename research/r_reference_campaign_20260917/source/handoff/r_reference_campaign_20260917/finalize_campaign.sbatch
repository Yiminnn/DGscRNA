#!/bin/bash
#SBATCH --job-name=Rref_finalize_delivery
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/Rref_finalize_delivery_%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PATH=/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin:$PATH
cd /fs/scratch/PCON0080/yimin/dgscrna
python handoff/r_reference_campaign_20260917/audit_design.py
python handoff/r_reference_campaign_20260917/audit_ptc_marker_seed_coverage.py --require-complete
python handoff/r_reference_campaign_20260917/plot_abstention.py
python handoff/r_reference_campaign_20260917/aggregate.py --require-complete
python -u handoff/r_reference_campaign_20260917/report_and_notebook.py
python -u handoff/r_reference_campaign_20260917/stage_delivery.py
rclone copy results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916/GBM_PTC_results_20260916 onedrive:work_od/share/dgscrna_GSE274546_TKU3186/01_report/GBM_PTC_results_20260916 --files-from results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916/R_reference_campaign_20260917_upload_files.txt --transfers 4 --checkers 8 --log-file results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916/R_reference_campaign_20260917_upload.log --log-level INFO
rclone check results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916/GBM_PTC_results_20260916 onedrive:work_od/share/dgscrna_GSE274546_TKU3186/01_report/GBM_PTC_results_20260916 --files-from results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916/R_reference_campaign_20260917_upload_files.txt --download --one-way --checkers 4 --log-file results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916/R_reference_campaign_20260917_remote_check.log --log-level INFO
python handoff/r_reference_campaign_20260917/record_delivery.py
rclone copyto results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/summary/DELIVERY_RECEIPT.json onedrive:work_od/share/dgscrna_GSE274546_TKU3186/01_report/GBM_PTC_results_20260916/results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/summary/DELIVERY_RECEIPT.json
