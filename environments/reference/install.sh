#!/usr/bin/env bash
# Usage: bash install.sh /path/to/conda-base/bin/python /path/to/new-reference-env
# Optional: DGSCRNA_DBSCAN_SOURCE=/path/to/dbscan_1.2.6.tar.gz (for offline CRAN access).
set -euo pipefail
if [[ "$#" != 2 ]]; then
  echo 'Usage: bash install.sh CONDA_BASE_PYTHON NEW_PREFIX_ROOT' >&2
  exit 2
fi
conda_python="$1"
reference_prefix="$2"
recipe_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
  echo 'These explicit locks support Linux x86-64 only.' >&2
  exit 2
fi
mkdir -p -- "$reference_prefix"
reference_prefix="$(cd -- "$reference_prefix" && pwd)"
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME R_HOME R_LIBS R_LIBS_USER R_LIBS_SITE LD_LIBRARY_PATH
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER='' R_LIBS_SITE=''
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$conda_python" -s "$recipe_dir/create_isolated.py" "$reference_prefix"
export PATH="$reference_prefix/r/bin:$PATH"
source_args=()
if [[ -n "${DGSCRNA_DBSCAN_SOURCE:-}" ]]; then
  source_args=(--source "$DGSCRNA_DBSCAN_SOURCE")
fi
"$reference_prefix/python/bin/python" -s "$recipe_dir/fetch_dbscan.py" \
  "$reference_prefix/sources/dbscan_1.2.6.tar.gz" "${source_args[@]}"
"$reference_prefix/r/bin/Rscript" --vanilla "$recipe_dir/install_dbscan_index64.R" \
  "$reference_prefix/sources/dbscan_1.2.6.tar.gz" "$reference_prefix/dbscan-source"
"$reference_prefix/r/bin/Rscript" --vanilla "$recipe_dir/check_isolation.R" \
  "$reference_prefix/r" "$reference_prefix/R_isolation.json"
"$reference_prefix/python/bin/python" -s "$recipe_dir/check_isolation.py" \
  "$reference_prefix/python" "$reference_prefix/Python_isolation.json"
"$reference_prefix/python/bin/python" -s "$recipe_dir/verify_downloads.py" \
  "$reference_prefix" "$reference_prefix/package_verification.json"
echo "Reference runtimes ready in $reference_prefix"
echo 'Install the matching DG-scRNA release wheel into the Python prefix with pip --no-deps.'
