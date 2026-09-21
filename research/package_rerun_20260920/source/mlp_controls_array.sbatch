#!/bin/bash
# Root controls campaign scheduling; this worker never submits another job.
#SBATCH --job-name=pkg_GBM_MLP
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/pkg_GBM_MLP_%A_%a.out
set -euo pipefail
task_gate=$1
task_here=/fs/scratch/PCON0080/yimin/dgscrna/handoff/package_release_20260920
task_python=$(/usr/bin/python3 -c 'import json,sys; g=json.load(open(sys.argv[1])); assert g["status"]=="release_and_pilots_verified"; print(g["runtime"]["python"])' "$task_gate")
task_output=$(/usr/bin/python3 -c 'import json,sys; from pathlib import Path; g=json.load(open(sys.argv[1])); t=json.load(open(sys.argv[2]))[int(sys.argv[3])]; print(Path(g["output_root"])/"extensions/MLP_controls"/t["sample"]/t["budget"]/t["control"])' "$task_gate" "$task_here/mlp_controls_tasks.json" "$SLURM_ARRAY_TASK_ID")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONNOUSERSITE=1
unset PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE DGSCRNA_REFERENCE_OUT DGSCRNA_EXAMPLE_OUT DGSCRNA_DL_CACHE_ROOT
cd /tmp
"$task_python" -s "$task_here/mlp_controls_adapter.py" --gate "$task_gate" --index "$SLURM_ARRAY_TASK_ID" --out "$task_output"
