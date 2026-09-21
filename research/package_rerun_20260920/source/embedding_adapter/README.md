# A1 adapter for the installed reference package

This research adapter consumes a **newly fitted, independently accepted packaged GBM preparation**, preserves the frozen reviewer A1 definitions, and invokes the installed package's original terminal DL backend. The public package and ongoing original A1 campaign are unchanged.

## Scientific contract

- 121 samples × HVG2000/HVG5000 × seven spaces = 1,694 representation units.
- Spaces: noDR, PCA2, FA2, ICA2 adaptive, Isomap2, direct-HVG UMAP2, TSNE2.
- Each space has KMeans/GMM at K=5/10/15/20/30/40 plus native-R HDBSCAN minPts50: 13 partitions, 22,022 conditions overall.
- The shared geometry is the exact float64 native-R scaled selected-HVG matrix. Cell/gene order is retained.
- **A1 UMAP uses scaled HVG directly → uwot**, cosine / n_neighbors30 / min_dist0.3 / seed42. This is a separate comparator from the package default **PCA30 → UMAP**.
- Scoring uses all retained RNA genes, original `wilcox_limma` DEG and density arithmetic, fixed CM2_glioma_other / mean cutoff. HDBSCAN noise0 remains a scored cluster.
- DL uses normalized selected-HVG RNA and the installed original architecture/training. Initial calls and terminal070/090 are retained; Lfine reporting uses terminal endpoints.
- ICA policy is copied byte-for-byte: parallel5000 → parallel50000 → parallel100000 → deflation100000, with the exact convergence checks. No truth or downstream score chooses a solver. Canonical failure records remain explicit.

`SOURCE_PROTOCOL_MANIFEST.json` records original hashes, unchanged Python functions and every allowed R substitution. `run.py` executes only the named numerical AST functions from the immutable copied source; it does not execute that source's old input/output globals. R modifications reverse exactly to the old bytes: input/marker paths, library isolation and lossless cell-ID parsing. `adaptive_ica.py` is unchanged.

## Runtime and inputs

`RUNTIME.json` records a dedicated installed-wheel venv. It inherits the original A1 numerical dependency environment read-only (including sklearn1.9.0); the minimal frozen package runtime_v3 is untouched. The wheel implementation hashes must match the accepted input run. All R calls use the pinned isolated package interpreter, `--vanilla`, and only `.Library`.

Required input: complete package run, matching independent `packaged_parity.json`, intact preparation/expression/DL hashes and the original marker source hash. Before each DL call, the scored cell file must match the accepted preparation cell-file hash. New caches stay inside each new representation directory; no old models, old representations, old DEG tables or old label arrays are loaded for fitting.

All execution requires SLURM, four DL threads, and OMP/BLAS/NUMBA thread limits1. Example within a four-CPU allocation:

```bash
"$A1_INSTALLED_PYTHON" -s handoff/package_release_20260920/embedding_adapter/run.py \
  --package-run "$ACCEPTED_PACKAGE_RUN" \
  --out "$NEW_PACKAGE_CAMPAIGN/embedding_adapter_v1/pilot" \
  --rscript "$PINNED_PACKAGE_RSCRIPT" \
  --spaces PCA2
```

Fresh execution is the default; existing outputs require explicit `--resume` and identical source/input/runtime configuration. Each space has an exclusive directory lock. Do not automatically remove a stale lock without checking its scheduler job.

## Validation scope

The bounded pilot refits TKU4163/HVG2000 for all seven spaces and 91 partitions from the new accepted full packaged pilot. `verify.py` independently compares geometry bytes, coordinates, partitions, initial calls, R DEG/density objects, terminal arrays, training splits/history and model weights against existing A1 results. It then invokes the separately frozen `evaluate_extension_lfine.py`, yielding 182 terminal-threshold Lfine rows.

Passed validation output:

`results/hvg_ptc_20260916_v1/package_reference_rerun_20260920/embedding_adapter_v1/verification_TKU4163_hvg2000/manifest.json`

The receipt reports `passed_exact`: seven representations, 91 terminal conditions and 182 Lfine threshold rows. Its SHA256 is `d80078665e498bddf4666fced6d3aa99d7a667824f7a56ee6f50db1163c0b07c`. This small pilot exercises ICA's successful parallel5000 branch; the adaptive fallback source is preserved but this is not a new fallback-run validation.

## Full-campaign work still required

No 1,694-unit submission is included. Full execution needs the root release gate, accepted new GBM preparations, a single geometry exporter per sample/budget before parallel representation jobs, pinned adapter/runtime/verification hashes, array scheduling/OOM monitoring, and independent acceptance. Existing original A1 jobs are neither adopted as new fits nor canceled. Patient-fold K selection and full-cohort comparisons must follow the existing Lfine protocol after all new candidates are complete; they are not established by this pilot.
