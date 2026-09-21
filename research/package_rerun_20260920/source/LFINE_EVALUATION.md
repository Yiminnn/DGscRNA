# Packaged GBM evaluation

The evaluator preserves the frozen v5 **compatible-target-set Lfine macro-F1**. Broad native marker labels map through the existing semantic parents into compatible Lfine classes. This is not strict one-to-one fine-type accuracy. The frozen native-R mapping includes the three previously documented Non-neuron corrections; no mapping, marker or model is selected using new outcomes.

- Full core: 726 sample/budget units × 4 routes × 16 libraries × 3 cutoffs = 139,392 terminal conditions. Each yields threshold 0.70 and 0.90, hence **278,784 evaluation rows**. Thresholds share a trained model.
- Class averaging: sample-present Lfine classes with support ≥20, excluding `Other`/`nan`; all cells remain in TP/FP/FN, exactly as v5.
- Cohort: fixed primary97 / sensitivity24 split. Primary97, all121, sensitivity24 and TKU3186 summaries retain expected/valid/unavailable counts. Patient means give equal patient weights after within-patient sample averaging.
- Terminal execution: trained, trained with one known class, no-op with all initial labels known, archived all-Undecided endpoint with no known labels, and structural insufficient split remain distinguished. They are accepted only with valid checked terminal receipts. Missing/corrupt/invalid endpoints are NA and never replaced by marker-only output.
- `BrainAtlas112`, `CARE_TME` and `UNION_all` carry the frozen author-reference overlap flag. The eight newer CellMarker libraries retain unspecified overlap as unknown; no independence is inferred.
- Output contains no coarse-label performance metric. Mapping audit source columns are ontology provenance only.

## Per-unit integration

Run after package fitting and parity verification, **inside the same SLURM allocation**:

```bash
"$RUNTIME_PYTHON" -s handoff/package_release_20260920/evaluate_core_lfine.py unit \
  --task-index "$SLURM_ARRAY_TASK_ID" \
  --campaign-root results/hvg_ptc_20260916_v1/package_reference_rerun_20260920
```

The task index refers to `rerun_inventory.core_tasks.json`. Reads `core/<sample>/<budget>/GBM/<sample>/<budget>/...`; writes `evaluation/units/<sample>/<budget>/metrics.csv.gz`, `manifest.json` and `COMPLETE`. No expression matrices or models are loaded.

Acceptance requires:

1. `COMPLETE` equals SHA256(`manifest.json`).
2. Manifest `status == "completed"`, `n_conditions == 384`, `n_valid == 384`, and task index/sample/budget match.
3. `outputs["metrics.csv.gz"]` matches the file hash.
4. `script_sha256` matches the pinned evaluator, and `frozen_semantics` matches `evaluate_core_lfine_sources.json`.

Pin both evaluator and source-lock files in the release gate. The source lock pins all nine provider/mapping/cohort/task files, verified before and after scoring. Invalid/unavailable unit evaluation writes its auditable receipt and returns exit code 2. A completed terminal with no scoreable Lfine class remains a valid terminal; its F1 is NA and each metric's available-sample count is explicit.

## Aggregation

Within SLURM:

```bash
"$RUNTIME_PYTHON" -s handoff/package_release_20260920/evaluate_core_lfine.py aggregate
```

Writes full metrics and summaries under `evaluation/`. Every expected condition is retained even when its unit has not been evaluated. `manifest.json` reports `partial` until all endpoints are valid. `metrics_hvg24` and `summary_hvg24` retain the fixed CM2_glioma_other / mean-cutoff comparison; full-gene marker summaries retain mean cutoff / UMAP2_HDBSCAN_R. Both thresholds are explicit; later notebook presentation should continue the established 0.90 main endpoint rather than select whichever wins.

The aggregate validates each per-unit receipt against the current frozen evaluator. Re-evaluate units after evaluator code changes. The `evaluation/validation/` subtree contains deliberate missing/tampered fixtures and a one-unit aggregation test; it is excluded from actual `evaluation/units/` aggregation.

## Validation

`verify_core_lfine.py` checks the full TKU4163 pilot (384 rows), NL022 all-gene pilot (2 rows), exact old terminal/metric agreement, five shared compact-notebook rows, 24 historical v5 endpoints, nontraining terminal status handling, missing and corrupted predictions, and a complete 278,784-row inventory with only one valid unit. Results and frozen hashes are in `evaluation/validation/verification.json`. No notebook or archived result was changed.
