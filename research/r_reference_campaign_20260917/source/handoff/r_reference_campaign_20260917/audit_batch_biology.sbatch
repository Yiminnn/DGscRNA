#!/bin/bash
#SBATCH --job-name=Rref_batch_biology
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/Rref_batch_biology_%A_%a.out
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /fs/scratch/PCON0080/yimin/dgscrna
task_code_dir=results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/execution_sources/${SLURM_JOB_ID}
mkdir -p "$task_code_dir"
cp handoff/r_reference_campaign_20260917/audit_batch_biology.py "$task_code_dir/audit_batch_biology.py"
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -u "$task_code_dir/audit_batch_biology.py" "$@"
