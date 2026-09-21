# DG-scRNA reference workflow

DG-scRNA 2.0.0rc1 packages the recovered **R preprocessing, clustering and marker scoring with the original Python DL/refinement**, behind one Python command. The numerical workflow is shared; the legacy simplified Python implementation remains separately available.

## Install

This is a GitHub prerelease on `align-r-reference`. Download the matching source archive and wheel from this repository’s Releases page. The verified runtime is Linux x86-64, Python 3.11, with pinned R/Python environments. `pip` installs the interface and backend files; it does not install R.

From the extracted source archive, install the fixed runtimes in a fresh directory:

```bash
bash environments/reference/install.sh /path/to/conda/bin/python /path/to/reference-env
/path/to/reference-env/python/bin/python -m pip install --no-deps /path/to/dgscrna-2.0.0rc1-py3-none-any.whl
/path/to/reference-env/python/bin/dgscrna doctor --rscript /path/to/reference-env/r/bin/Rscript
```

On an HPC system, run installation and computation in an allocation. Portable SLURM examples are in `examples/reference/`.

## Run

Provide an **already-QC**, single-sample 10x count matrix directory and a marker library. The GBM preset does not repeat raw-droplet QC, DoubletFinder or cross-patient CCA.

```bash
dgscrna run --counts sample_10x/ --markers markers.tsv \
  --preset gbm-reference --sample sample1 --out results/sample1 \
  --rscript /path/to/reference-env/r/bin/Rscript
```

The default is VST2000 → PCA30 → uwot UMAP2 → R HDBSCAN → Seurat RNA DEGs → marker density with mean cutoff → original DL/refinement. Geometry and DL use the selected genes; RNA DEG scoring uses all eligible genes. `--features all` uses all genes detected in at least three cells. Other supported budgets are 500/1000/3000/5000. `--route all --cutoff all` runs the four original routes and three scoring cutoffs. This preset is a reference configuration, not a claim of universal optimality.

Markers can be nested JSON `{library: {cell_type: [genes]}}` or a long TSV/CSV:

```text
library	cell_type	gene
brain	Astrocyte	GFAP
brain	Astrocyte	AQP4
brain	Oligodendrocyte	MBP
```

The run returns `annotations.csv.gz` (initial/final labels, confidence, cluster and DL status for each condition), `embedding.csv`, `condition_summary.csv`, and `run_manifest.json`. Intermediate artifacts and logs remain available. Multiple libraries produce separate conditions; the tool does not use reference truth to pick a winning annotation. Valid no-op outcomes are distinguished from actual training.

Python uses the same runner:

```python
from dgscrna import run_reference
result = run_reference(counts="sample_10x", markers="markers.tsv",
                       out="results/sample1", sample="sample1",
                       rscript="/path/to/reference-env/r/bin/Rscript")
```

See [the reference guide](docs/reference_workflow.md) for h5ad inputs, resume, staged HPC execution, provenance and limitations. Prepared-input PTC CCA code is included as an advanced backend; raw PTC preprocessing and end-to-end PTC portability are not certified by the GBM fixture.

## Reproducibility and legacy code

`research/reference_examples/` preserves the independently verified TKU4163 fresh-fit fixture. The package carries backend source hashes and tests that preserve the original numerical bodies; release-specific fresh-fit checks are reported in release notes. Fitting never requires author reference labels.

Existing `dgscrna.core` functions remain the **legacy simplified Python implementation**, available with `pip install 'dgscrna[legacy]'`. They are not interchangeable with the R reference and are not used by the new `dgscrna run` command.

GPL-3.0. Report issues at https://github.com/Yiminnn/DGscRNA/issues.
