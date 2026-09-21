# Historical Table 2 saved-label reconciliation

These unchanged scripts document the previously verified historical replay. They are **site-path source snapshots**, not a portable Table 2 refit. No PTC numerical job has been run for this example package; a fresh replay remains gated on completion of GBM.

Original execution order is `check_paper_baseline.py` → `reconcile_historical_metrics.py` → `verify_metric_R.R` → `check_accuracy_consistency.py`. All require SLURM. The first script also audits supplementary-label/branch provenance and therefore has more data dependencies than a single metric calculation. Its literal paths must be redirected to a fresh output directory before any future execution; running this unchanged snapshot would target the old baseline directory and is not recommended.

Required historical inputs:

- `tcr/rawdata/data_with_validation.csv`: original S3 final native/detailed labels and original validation flags.
- `tcr/rawdata/data_with_validation+3_cell_types 2.csv` and `metadata.txt`: original branch provenance and sample mapping.
- `tcr/ptc_val/scripts/DGscRNA-Vignette.rmd`: original label-collapse definitions; the historical broad T list includes NK labels.
- `tcr/ptc_val/scripts/DGscRNA-Share/R/source.R` and `ptc_recovery/inventory_workspace/object_01.meta.csv`: branch auditing.
- Original supplementary S3 workbook and manuscript V16 DOCX under `paper/submission_v16/`.
- Archived `MLmetrics_1.1.1_Classification.R` included here: verifies that the original default F1 used class0 as positive.

The already extracted archive is under `results/hvg_ptc_20260916_v1/ptc_recovery/archive/`; its original source is the user archive `work/_archives/tcr.tar.gz`. Identify final endpoints by sample/barcode and content, not changed filenames. S2 (`tcr/scripts/annotations.xlsx`) is a separate endpoint from S3 and must remain separate.

Existing verified outcome: **31/36 numeric manuscript entries agree at reported precision; all eight DG-scRNA F1/AUC entries agree.** This is not full Table 2 reproduction. SCINA Primary tumor F1 rounds to0.9366 from stored labels rather than0.9367. Four methods' published overall Accuracy entries do not equal ordinary binary accuracy for the recovered endpoint; their source remains unrecovered. Do not change labels or select a different endpoint to force agreement. Historical SignacX predicts0 T cells and its literal manuscript `None` entries stay unchanged.

The supporting site records are `ptc_paper_baseline/historical_metric_reconciliation_manifest.json`, `paper_Table2_historical_source_comparison.csv`, and `paper_accuracy_consistency_manifest.json`. Replaying saved final labels proves this endpoint accounting; it does not recover missing historical DL initialization/weights or establish exact fresh-model reproduction.
