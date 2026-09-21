#!/usr/bin/env bash
#SBATCH --job-name=dgscrna-gbm
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=dgscrna-gbm-%j.log
# Supply DGSCRNA_ENV, DGSCRNA_COUNTS, DGSCRNA_MARKERS, DGSCRNA_OUT and DGSCRNA_SAMPLE.
# Add site-specific account/partition flags to sbatch. Increase memory for large samples.
set -euo pipefail
: "${DGSCRNA_ENV:?Set the installed reference runtime root}"
: "${DGSCRNA_COUNTS:?Set the 10x input directory}"
: "${DGSCRNA_MARKERS:?Set the JSON or long-TSV marker file}"
: "${DGSCRNA_OUT:?Set a new result directory}"
: "${DGSCRNA_SAMPLE:?Set a unique sample name}"
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH PYTHONHOME R_HOME R_LIBS R_LIBS_USER R_LIBS_SITE LD_LIBRARY_PATH
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER='' R_LIBS_SITE=''
export PATH="$DGSCRNA_ENV/r/bin:$PATH"
"$DGSCRNA_ENV/python/bin/dgscrna" run \
  --counts "$DGSCRNA_COUNTS" --markers "$DGSCRNA_MARKERS" \
  --preset gbm-reference --sample "$DGSCRNA_SAMPLE" \
  --features "${DGSCRNA_FEATURES:-2000}" --route "${DGSCRNA_ROUTE:-UMAP2_HDBSCAN_R}" \
  --library "${DGSCRNA_LIBRARY:-all}" --cutoff "${DGSCRNA_CUTOFF:-mean}" \
  --rscript "$DGSCRNA_ENV/r/bin/Rscript" --out "$DGSCRNA_OUT" --require-slurm
