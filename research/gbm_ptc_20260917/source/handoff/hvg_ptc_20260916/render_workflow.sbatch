#!/bin/bash
#SBATCH --job-name=gbm_hvg_workflow
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/gbm_hvg_workflow_%j.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
OUT=results/hvg_ptc_20260916_v1/figures
mkdir -p "$OUT"
dot -Tsvg handoff/hvg_ptc_20260916/workflow.dot -o "$OUT/decision_tree.svg"
dot -Tpdf handoff/hvg_ptc_20260916/workflow.dot -o "$OUT/decision_tree.pdf"
dot -Tpng -Gdpi=140 handoff/hvg_ptc_20260916/workflow.dot -o "$OUT/decision_tree.png"
