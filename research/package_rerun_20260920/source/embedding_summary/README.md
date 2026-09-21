# Packaged A1 patient-fold comparison

This summary runs only after the new complete extension campaign (2,617 accepted tasks, including all 1,694 A1 units) and all 726 core runs are accepted. It consumes the new shared Lfine metrics and same-budget new core anchors. It does not fit, reuse historical selected K, or substitute partition scores for terminal annotation.

The numerical selection and independent verification bodies come from the frozen `handoff/reviewer_completion_20260920/embedding_lfine_v1/select_lfine.py` and `verify_selection.py`. `SOURCE_PROTOCOL_MANIFEST.json` records their exact hashes and the complete substitutions. Only endpoint-count assertions change from three stages to the two terminal thresholds; marker-only diagnostics are omitted. No arithmetic, ranking, inference or missingness rule changes.

- Fixed author patient folds; primary97 and all121 cohorts evaluated separately.
- Within-patient sample mean, then equal training-patient mean terminal090 Lfine.
- KMeans/GMM K=5/10/15/20/30/40; exact ties choose smaller K. HDBSCAN remains fixed.
- The held-out fold does not choose K. Terminal070 uses the same terminal090-selected K.
- All 21 space/clusterer contrasts remain; 10,000 paired patient bootstrap draws, seed20260920; Wilcoxon and Holm21 per cohort/budget.
- Fixed same-budget original PCA30→UMAP→HDBSCAN anchor. A1 UMAP instead starts directly from scaled HVGs.
- Truth-only scored-class eligibility must agree across every candidate and anchor. Zero-class entries remain explicit; they are not missing fits.

Expected outputs: 44,044 candidate terminal rows, 484 anchor rows, 420 fold choices and 84 paired contrasts. The independently implemented frozen verifier must reproduce every selection, denominator and statistic before completion is written.

```bash
"$A1_INSTALLED_PYTHON" -s handoff/package_release_20260920/embedding_summary/run.py \
  --campaign "$UNIFIED_EXTENSION_GATE" \
  --aggregate "$NEW_PACKAGE_CAMPAIGN/evaluation/extensions" \
  --update-notebook
```

Run through SLURM. Without `--update-notebook`, only a local notebook candidate is written after the complete selection passes. With it, a shared lock and backup protect the existing notebook; only its DGCyTOF subsection changes. Core, reviewer and PTC sections remain intact. No HTML or remote synchronization occurs.

Validation: incomplete-campaign rejection and two full-roster synthetic regressions passed in SLURM. Both regressions cover 44,044 terminal candidate rows, 484 anchors, 420 fold choices, 84 contrasts, exact K ties, held-out-label perturbation, a truth-ineligible patient, parity to the original three-stage producer's terminal outputs, and notebook candidate preservation. Real campaign aggregation remains pending actual completion.

The quantized regression preserved a numerical boundary in the original independent verifier: its `groupby.mean` and the producer's `Series.mean` differ near machine precision, changing Wilcoxon tie ranks (p difference 5.746e-12, Holm difference 3.984e-11). `canonical_inference.py` independently verifies all patient means and non-p statistics at the unchanged tolerance, rebuilds canonical candidate means exactly from selected sample rows, and then independently recomputes Wilcoxon/Holm from losslessly saved patient vectors. Both quantized and continuous regressions reproduce the producer's p values exactly. No producer arithmetic, ranks or tolerance changes. The original failing fixture and log remain preserved.

Receipts: `results/hvg_ptc_20260916_v1/package_reference_rerun_20260920/embedding_summary_validation/{quantized_v3,continuous_v3}/verification.json`. `independent_review.json` binds the reviewed source and both regressions. None is a completion claim for real A1 fits or cohort results.
