#!/bin/bash
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=7-00:00:00
#SBATCH --job-name=reviewA_dispatch
set -euo pipefail
cd /fs/scratch/PCON0080/yimin/dgscrna
exec /users/PCON0080/yimin/miniforge3/bin/python "$1"
