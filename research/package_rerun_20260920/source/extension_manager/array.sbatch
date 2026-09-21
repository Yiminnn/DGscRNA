#!/bin/bash
# Root activates only a reviewed extension campaign gate.
#SBATCH --job-name=pkg_GBM_extensions
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --no-requeue
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/pkg_GBM_extensions_%A_%a.out
set -euo pipefail
campaign_gate=$1
task_python=$(/usr/bin/python3 -c 'import json,sys; c=json.load(open(sys.argv[1])); assert c["status"]=="reviewed_extension_campaign"; r=json.load(open(c["release_gate"]["path"])); print(r["runtime"]["python"])' "$campaign_gate")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONNOUSERSITE=1
unset PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE
cd /tmp
"$task_python" -s /fs/scratch/PCON0080/yimin/dgscrna/handoff/package_release_20260920/extension_manager/worker.py --campaign "$campaign_gate"
