#!/bin/bash
# Root submits only after a published release and four accepted pilots.
# Usage: sbatch --array=0-725%N core_array.sbatch /absolute/RELEASE_GATE.json
#SBATCH --job-name=pkg_GBM_core
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/pkg_GBM_core_%A_%a.out
set -euo pipefail
task_gate=$1
task_python=$(/usr/bin/python3 -c 'import json,sys; x=json.load(open(sys.argv[1])); assert x.get("status")=="release_and_pilots_verified", "release gate pending"; print(x["runtime"]["python"])' "$task_gate")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONNOUSERSITE=1
unset PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE
cd /tmp
"$task_python" -s /fs/scratch/PCON0080/yimin/dgscrna/handoff/package_release_20260920/core_task.py --gate "$task_gate"
