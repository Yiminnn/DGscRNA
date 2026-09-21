# Validated native R GBM reference example

This separate research example runs **raw counts → native R preprocessing → UMAP2/HDBSCAN → original density marker scoring → terminal DL/refinement**. It does not call the simplified Python package API. The fixed TKU4163 example uses 2,000 HVGs and genuinely trains the final DL stage; it is a reproducibility fixture, not a claim that 2,000 HVGs or this route is optimal.

`SOURCE_DERIVATION.json` records every original source hash and every IO-only adaptation. The R numerical bodies are unchanged; `terminal.py` and `legacy_refine.py` are byte-identical to the frozen reference. The original all-RNA DEG scorer still uses all retained genes; HVG selection controls geometry and the normalized DL input. Both terminal confidence endpoints, 0.90 and 0.70, are retained. Reference labels are never opened by the fitter.

## Data and execution

The three input files are the author-retained TKU4163 counts from GSE274546: `matrix.mtx` (genes × cells), `genes.tsv`, and `barcodes.tsv`. Their exact hashes and expected dimensions are in `fixtures/TKU4163.json`. The full frozen `markers/libraries.json` is required with its recorded SHA256; only the predeclared CM2_glioma_other/mean arm is sent to DL. Raw count and marker files are external data dependencies and are not included here.

Run within a SLURM allocation with 4 CPUs and 32 GB RAM:

```bash
python run_gbm.py --counts-dir /data/TKU4163 \
  --markers /data/markers/libraries.json \
  --rscript /env/reference-r/bin/Rscript --output /results/fresh-example
```

The output directory must not exist. No earlier predictions, normalized matrices, R objects or DL cache are copied. Input conversion must exactly match the frozen binary counts/cell/gene hashes before R starts. Two R DEG workers and fixed seed42 are used. The final predictions are under `GBM/TKU4163/hvg2000/UMAP2_HDBSCAN_R/terminal/L00_mean/`.

To compare an independently fitted example with an existing frozen reference:

```bash
python verify_gbm.py --output /results/fresh-example \
  --reference /reference/GBM/TKU4163/hvg2000
```

The runner automatically checks the frozen semantic fingerprints in `fixtures/TKU4163_expected_outputs.json` after fitting; no prior result directory is needed for that check. The optional explicit-reference verifier additionally checks feature order, normalized DL matrix, embeddings, clusters, all initial marker arms, DEG/density R objects, DL probabilities, train/validation cell indices, terminal labels and history. Serialized R-object byte hashes may contain runtime metadata and are not used as a substitute for numerical comparisons.

## Isolated environment recipe

`environment/*-linux-64.explicit.txt` pins the installed package dependency closures for two fresh Conda prefixes; this is Linux x86-64 specific. Use a fresh package cache and copy installation to avoid previously modified hardlinked package-cache files. Installed environments and result caches are separate. The R prefix additionally needs the source-built dbscan1.2.6.9001 integer-index fix. `install_dbscan_index64.R` accepts the verified upstream dbscan1.2.6 source archive and a fresh extraction directory. This restores the reference fix used for large PTC matrices; no PTC calculation is part of this GBM example.

```bash
mkdir /env/reference
/conda-base/bin/python environment/create_isolated.py /env/reference
export PATH=/env/reference/r/bin:$PATH
unset R_LIBS R_LIBS_USER R_LIBS_SITE PYTHONPATH
/env/reference/r/bin/Rscript --vanilla environment/install_dbscan_index64.R \
  /data/dbscan_1.2.6.tar.gz /results/dbscan-source
```

The installer uses Conda’s documented `search_path=()` API override and asserts that the effective cache list contains only the fresh cache. This avoids user/site configuration merging an old cache into an intended isolated installation. It does not change global configuration.

Set `PYTHONNOUSERSITE=1` and use the new Python/R executables. The isolation checks record actual loaded numerical-package paths, so a nominal new prefix cannot silently reuse the old R/Python libraries. The fixed GBM example passed an actual clean-prefix installation and fresh fit in SLURM7447099. All loaded R/Python numerical libraries resolved inside the new prefixes, and331 Conda package archive/build records passed independent checksum checks. Fresh outputs exactly matched the historical reference, including DEG/density R objects and actual DL probabilities, split, history and both final-label endpoints. This certifies this Linux x86-64 fixture and recorded dependencies; other datasets, hardware and environments have not been certified by this smoke test.

The current site data are under `data_bench/GSE274546/mtx/TKU4163/`; frozen markers and reference results are under `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/`. New job receipts and parity proofs are under `results/hvg_ptc_20260916_v1/reviewer_completion_20260920/reference_examples/` at the workspace root. These relative site paths are provenance, not a portable download service.

A separately prepared 2.25 MB count-and-marker fixture archive is available in the site result directory `reference_examples/fixture_bundle_v1/`, with member hashes and source provenance. It contains no final annotations, model weights or fitted cache.

Validation receipts: `current_env_7447043/parity_v2.json`, `clean_env_7447099/parity.json`, `isolated_7447099/{R,Python}_isolation.json`, and `pristine_package_verification_7447099.json` in the site result directory. Attempts7447052/7447064 stopped during installation because the old shared cache was unexpectedly selected; those failed attempts are preserved and the successful installation explicitly excluded all user/site cache configuration. No global cache or environment was repaired in place.
