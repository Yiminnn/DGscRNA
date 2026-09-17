#!/bin/bash
#SBATCH --job-name=ptc_harmony_isolation
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/ptc_harmony_isolation_%j.out
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -m pip install --no-deps --target results/hvg_ptc_20260916_v1/ptc_recovery/vendor_harmony harmonypy==0.0.10
PYTHONPATH=/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/ptc_recovery/vendor_harmony /fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python -c 'import harmonypy; print(harmonypy.__version__, harmonypy.__file__); assert harmonypy.__version__ == "0.0.10"'
