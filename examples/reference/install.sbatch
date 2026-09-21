#!/usr/bin/env bash
#SBATCH --job-name=dgscrna-install
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=dgscrna-install-%j.log
# Submit from the repository root, adding your scheduler's account/partition flags:
# sbatch --export=ALL,CONDA_BASE_PYTHON=/conda-base/bin/python,DGSCRNA_ENV=/env/new-reference \
#   examples/reference/install.sbatch
set -euo pipefail
: "${CONDA_BASE_PYTHON:?Set CONDA_BASE_PYTHON to Conda base's Python executable}"
: "${DGSCRNA_ENV:?Set DGSCRNA_ENV to a fresh writable runtime directory}"
cd -- "${SLURM_SUBMIT_DIR:?Submit this script with sbatch from the repository root}"
bash environments/reference/install.sh "$CONDA_BASE_PYTHON" "$DGSCRNA_ENV"
# Install the release wheel afterward, or pass its absolute path in DGSCRNA_WHEEL.
if [[ -n "${DGSCRNA_WHEEL:-}" ]]; then
  "$DGSCRNA_ENV/python/bin/python" -m pip install --no-deps "$DGSCRNA_WHEEL"
fi
