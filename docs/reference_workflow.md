# Reference workflow

DG-scRNA 2.0.0rc1 exposes the recovered R workflow and historical Python DL/refinement through one command. The `gbm-reference` preset uses Seurat for normalization, VST features, PCA, UMAP and marker scoring; R `dbscan` for HDBSCAN; and the recovered Python training protocol for terminal labels. The older simplified Python workflow remains a separate implementation.

## Install once

The reproducible runtime consists of **two Conda prefixes**, R and Python. Installing the Python wheel alone does not install R or Seurat. The explicit locks in `environments/reference/` support Linux x86-64 and preserve the dependency closure used by the clean-environment reference fixture.

From a release source checkout, run:

```bash
bash environments/reference/install.sh \
  /path/to/conda-base/bin/python /path/to/new-reference-env

/path/to/new-reference-env/python/bin/python -m pip install --no-deps \
  /path/to/dgscrna-2.0.0rc1-py3-none-any.whl
```

Use Conda base's Python for installation: it must provide the `conda` API. The installer uses a fresh package cache and copy installation, excludes user/site Conda configuration, verifies the package archives, and builds the recorded `dbscan` index-width fix from checked upstream source. It does not change global Conda settings. Installation needs network access and a writable directory with enough space for two runtimes and their package archives. An existing `dbscan` archive can be supplied with `DGSCRNA_DBSCAN_SOURCE`; its SHA256 must match.

On an HPC cluster, use `examples/reference/install.sbatch` and supply your account/partition via `sbatch`. The original fixture was independently installed and fitted on Linux x86-64; this does not certify other operating systems or all datasets. Runtime manifests record the actual versions used. These pins include the historical PyTorch build; the reference fit runs on CPU and does not require a GPU allocation.

## Run GBM

```bash
/path/to/new-reference-env/python/bin/dgscrna run \
  --counts /data/sample_10x \
  --markers /data/markers.json \
  --preset gbm-reference \
  --sample sample01 \
  --rscript /path/to/new-reference-env/r/bin/Rscript \
  --out /results/sample01
```

For SLURM, use `examples/reference/run_gbm.sbatch`. Each output directory belongs to one configuration. Use a different output directory when changing the feature budget or other parameters. `--resume` resumes the same configuration only after verifying the recorded stage artifacts; `--stop-after prepare` or `--stop-after score` produces an explicitly incomplete run. Only the default `--stop-after complete` attempts terminal refinement.

The Python interface uses the same backend:

```python
from dgscrna import run_reference

run_reference(
    counts="/data/sample_10x",
    markers="/data/markers.json",
    out="/results/sample01",
    sample="sample01",
    features="2000",
    route="UMAP2_HDBSCAN_R",
    library="all",
    cutoff="mean",
    rscript="/path/to/new-reference-env/r/bin/Rscript",
)
```

## Input contract

- **Counts:** nonnegative integer RNA counts from cells that have already passed the intended study QC. The 10x directory contains `matrix.mtx`, `barcodes.tsv` and `genes.tsv` or `features.tsv`; gzip variants are accepted. The matrix is genes × cells. Cell IDs and gene names must be unique and genes must match the marker namespace. Resolve duplicate symbols explicitly before fitting.
- **AnnData:** `.h5ad` input requires `--counts-layer NAME`, or explicit `--counts-layer X` if `.X` contains raw counts. A normalized expression matrix is not a count matrix. AnnData support requires the optional `anndata` dependency.
- **Markers:** JSON with the structure below, or a long TSV with columns `library`, `cell_type`, `gene`. All libraries run by default; `--library NAME` selects one. Library names and cell-type names are preserved. The shown genes demonstrate syntax only; use a complete, documented marker panel for analysis.
- **Size:** the reference PCA30 workflow requires more than 31 retained cells and at least 31 usable features. Genes detected in fewer than three retained cells are filtered.
- **Truth labels:** not an input to fitting. Evaluate saved final predictions separately using a prespecified label mapping and denominator. A marker vocabulary may be coarser than the author's labels.

```json
{
  "brain_context": {
    "Astrocyte": ["AQP4", "GFAP"],
    "Oligodendrocyte": ["MBP", "PLP1"]
  }
}
```

```text
library	cell_type	gene
brain_context	Astrocyte	AQP4
brain_context	Astrocyte	GFAP
brain_context	Oligodendrocyte	MBP
brain_context	Oligodendrocyte	PLP1
```

## Preset and comparison options

| Setting | `gbm-reference` default | Available comparison |
|---|---|---|
| Cell QC | Use the supplied retained cells | Full raw-cell QC is not implemented by this preset |
| Gene filter | Detected in at least 3 cells | Recorded in the input manifest |
| Normalization | Seurat LogNormalize, scale factor 10,000 | Fixed reference setting |
| Features | VST 2,000 | `--features 500`, `1000`, `2000`, `3000`, `5000`, `all` |
| Batch correction | Single sample; no correction | PTC integration is a separate workflow |
| Geometry | Scale → PCA30 → uwot UMAP2 | PCA30 or UMAP2 clustering routes |
| UMAP | Cosine, 30 neighbors, min_dist 0.3, seed 42 | Fixed reference setting |
| Clusterer | R HDBSCAN on UMAP2 | `PCA30_HDBSCAN_R`, `PCA30_SNN`, `UMAP2_SNN`, or `--route all` |
| Marker seed cutoff | `mean` | `--cutoff none`, `0.5`, or `all` |
| Terminal endpoint | Recorded 0.90 and 0.70 confidence endpoints | Preserve both; choose the endpoint before evaluation |

`--features all` means all genes retained after the detection filter, not the entire genome. It skips HVG selection but still uses PCA and UMAP; it is **not** a no-dimensionality-reduction experiment. HVG selection controls geometry and the normalized DL input. GBM marker scoring uses all retained normalized RNA genes, including genes outside the HVG set. The exact selected, geometry, scoring and DL gene lists are saved separately.

The 2,000-gene default is the recovered reference setting. Neither this default nor the packaged UMAP/HDBSCAN route is a claim of optimal performance. Comparisons must retain the same cells, markers, endpoint and evaluation mapping.

`--seed` controls R preparation, including PCA and UMAP. The historical DL model initialization and train/validation split retain seed 42; changing `--seed` is a geometry sensitivity experiment, not a joint R-and-DL random-seed experiment.

## Read the result

Start with these files in the output root:

| File | Contents |
|---|---|
| `annotations.csv.gz` | One row per cell × route × marker library × cutoff; includes `initial`, `final090`, `final070`, confidence, lineage, cluster and DL status |
| `embedding.csv` | Actual native-R UMAP coordinates from this run, indexed by cell ID; shared across marker arms |
| `condition_summary.csv` | One row per terminal arm, with status and prediction-file path/hash |
| `run_manifest.json` | Package/configuration provenance, terminal and trained-arm counts, and export hashes |
| `GBM/` | Full preparation, scoring and terminal artifacts |
| `logs/` | Per-stage execution logs |

Use terminal predictions, not initial marker calls, as the DG-scRNA annotation endpoint. Each arm also saves `predictions.csv.gz` and `terminal_manifest.json` under `GBM/<sample>/<condition>/<route>/terminal/<arm>/`. The manifest identifies its marker library and cutoff, input hashes, terminal validity, whether training actually ran, and any exact-input result reuse. Preparation and scoring manifests retain the gene lists, embeddings, clustering and marker provenance. Join coordinates by cell ID; a cell has repeated annotation rows when several conditions run.

| `dl_status` | Meaning |
|---|---|
| `trained` | DL trained and processed the unresolved pool |
| `trained_single_known_class` | DL ran with one known seed class; interpretation is limited |
| `no_op_all_initially_known` | No unresolved cells; terminal labels retain the initial calls |
| `no_known_labels_archived_Undecided_terminal` | No known seeds; cells remain unresolved and no model was fitted |
| `structural_insufficient_known_split` | Too few known cells for the historical split; unresolved cells receive `Unknown` |

An explicitly completed no-op is a terminal result with its status recorded. It is not evidence that a model trained. A failed or incomplete stage is not a valid annotation result. Known seed labels remain unchanged by the historical refinement; inspect coverage and unresolved cells as well as classification scores.

Fitting does not emit accuracy or select markers from author labels. Saved labels can be evaluated afterward; keep Lfine concordance, strict subtype classification and clustering metrics distinct. The GBM report uses the frozen Lfine mapping. PTC S2 and S3 are distinct final-label sources and must retain their respective definitions.

## Current scope

The general GBM entry point packages the already-QC, single-sample reference path. The fixed research example under `research/reference_examples/` remains the numerical parity fixture, with independent clean-environment verification. Full raw-cell QC, DoubletFinder and a raw-count PTC CCA input contract are not part of the GBM preset. Prepared PTC research objects are not a substitute for a validated public raw-data interface.

The requested PTC grouped refits use MT-1/MT-2/N-1/N-2 and TU-1/TU-2/T-1/T-2. These are distinct from the historical saved eight-sample CCA context. Original marker selections, label simplification and final-label endpoints need their own baseline verification before ablations. A separate PTC release will document that validated input contract rather than silently run the single-sample GBM preset on pooled cells.

For errors, retain the run directory, stage logs and manifests. Use a fresh output directory for a changed analysis. Existing results are never evidence of a successful new fit unless exact-input reuse is explicitly recorded.
