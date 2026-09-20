# Paper comparator audit for the reviewer completion plan

Audit date: 2026-09-20. Read-only inspection of existing scripts, paper extracts, and saved summary tables. No new experiment, matrix load, plotting, or result modification was performed.

## Confirmed decision

The supplied V16 paper reports **zero T-cell detections for SignacX, not scType**. This is supported independently by the paper's Table 2, the archived writing draft, Reviewer 3 B3/B4, and the S3 reconciliation tables. There is no evidence here permitting scType to be changed to zero.

The user subsequently confirmed: **“指 SignacX，按原稿保留 0Tcell”**. Preserve the original paper's SignacX zero-T result and retain the actual scType results. This is a settled decision, not a pending clarification. New audit/rerun results must be separate records and must not overwrite the original paper/S3 record or be forced to zero.

| Source | scType | SignacX |
|---|---|---|
| V16 Table 2, overall | F1 0.9392; AUC 0.9271; Accuracy 0.8823 | `None` in all contexts and metric rows |
| Archived writing draft | “scType (F-1: 93.92%)” | “SignacX did not predict any cell as T cell.” |
| Literal supplied S3 flags, all 92,404 cells | 37,468 predicted positive | 0 predicted positive |
| Literal S3 flags, NMT / TTU | 27,747 / 9,721 predicted positive | 0 / 0 predicted positive |
| Reviewer 3 | B4 names scType as strongest reported comparator | B3 explicitly challenges the zero-T result and asks to verify, correct, or remove |

“Zero” here means no cells called T under the archived output/flag, not zero accuracy, not zero AUC, and not proof the program failed to produce predictions. The paper printed `None`; keep that exact reported entry in the historical table, with `0 predicted T cells` in a descriptive field. Literal all-negative binary evaluation has T-positive F1 = 0 and AUC = 0.5; these are separately recomputed diagnostics, not the paper's reported entries.

## Why the saved scType metrics differ across tables

These are distinguishable endpoints, not alternative numbers from which to select the most favorable result:

1. `ptc_paper_baseline/paper_Table2_historical_source_comparison.csv` reproduces all eight reported scType F1/AUC values at four decimals. It uses the archived CSV's original `validation_t_cell` (35,727 positives, matched to any filtered TCR contig) and the original scType binary flag. The historical MLmetrics F1 default treats class **0** as positive. Overall reconstructed F1 is 0.9392454284 and binary-call AUC is 0.9271423511.
2. `ptc_paper_baseline/paper_S3_literal_metrics.csv` uses supplied S3 `T_cell` (36,184 positives) with positive-class-1 F1. Its overall scType T-positive F1 is 0.9679302667. It should not be described as reproducing the published F1 definition.
3. `paper_claim_validation_20260917/PTC_comparator_replay/all_assay_metrics.csv` normalizes archived native label spellings through a frozen ontology and evaluates **strict T**, separating NK/NKT from T, against several named TCR definitions. Productive-TCR strict-T overall scType F1 is 0.9126755796. It is a re-evaluation of cached predictions, not a new scType fit and not the historical broad flag.
4. The V16 scType Accuracy entry 0.8823 does **not** match ordinary binary accuracy from the recovered original endpoint (0.9266157309). Preserve it as paper-reported and mark its definition/source unresolved. Do not copy it into a new ordinary-accuracy table.

## Original implementation provenance

The original archived Rmd at `results/hvg_ptc_20260916_v1/ptc_recovery/archive/tcr/ptc_val/scripts/DGscRNA-Vignette.rmd`:

- Lines 691–702: SignacX receives the integrated object, with integrated `counts` assigned from integrated `data`; calls `SignacFast(..., do.normalize = T, graph.used = 'snn')`, then `GenerateLabels`.
- Lines 714–763: scType loads original functions and applies `sctype_score(scRNAseqData = integrated_data@assays$integrated@data, scaled = F, gs = markers, gs2 = NULL)`.
- Lines 795–812: scType sums scores within Seurat clusters, takes the top-scoring type, and assigns those labels.
- Line 820: archived SignacX annotations use `signacx.celltypes$CellStates`.

The recovered original S3 SignacX native labels include epithelial, endothelial, unclassified, DC, fibroblast, classical monocyte, and neutrophil calls, with no lymphoid calls. This is an actual populated prediction vector. Its failure cause is **unestablished**. Existing comments suggesting Assay5 or missing graph are hypotheses and must not be presented as a proven explanation.

Read-only comparison of the original Rmd and `handoff/ptc/run_signacx.R` establishes that the later nonzero output is a **different fit**: original pooled integrated expression plus `SignacFast`; later per-sample raw 10X counts, RNA normalization, 2,000 variable genes, PCA30/FindNeighbors, v3-assay compatibility option, and `Signac`. The old output also uses `CellStates`, so its zero cannot be attributed merely to evaluating `CellTypes` instead of `CellStates`. Several inputs and execution choices changed together; this does not identify which change caused the nonzero result. No claim that the historical failure was fixed by one particular change is justified by the existing evidence.

## Existing SignacX follow-up is separate evidence

Completed replay job 7398160 verified the cached SignacX 2.2.5 concatenated output against all eight per-sample files and joined 92,404 of the 107,545 source cells to the historical cohort. It made **zero new fits**. It separately analyzes explicit T-prefixed `CellStates` and coarse `CellTypes == TNK`; TNK includes NK and is not a strict-T endpoint.

The corrected cached strict-T results are nonzero (patient-mean F1 approximately 0.8959 in NMT and 0.8312 in TTU). They do not erase the original paper's zero-T result. They also cannot be suppressed while asserting Reviewer 3 B3 was resolved by retaining only the historical zero result. Historical inputs/pooling/reference information differ; the existing replay is not a fully matched-information benchmark.

## Proposed completion-plan treatment

1. Freeze a **paper-reported** comparator table using the V16 entries exactly, as the user confirmed: scType 0.9392 F1; SignacX `None`, 0 predicted T cells. Keep native S3 labels and literal flags immutable.
2. Include a **historical endpoint reproduction** table with the recovered metric definitions, source hashes, per-group cells, and remaining Accuracy discrepancy. This completes historical provenance, not fairness validation.
3. Include a **revision comparison** with a single declared cell universe, TCR definition, label mapping, Unknown handling, patient unit, and allowed marker/parameter-selection opportunities. Cached outputs can be reused only where their inputs support that protocol. Record unmatched information explicitly.
4. Address Reviewer 3 B3 with a provenance appendix: original input assay, pooled versus per-sample processing, native vocabulary, zero-T count, software/source availability, and the separate later result. Preserve the historical zero table exactly; do not use that table alone to claim a fair superiority comparison. An additional modern run is needed only if source/configuration verification shows the existing later output cannot support the declared protocol. Any such run is separately labeled and its actual result is retained. If a fair supported SignacX comparison is not feasible, exclude it from the headline ranking with an explanation while retaining historical provenance. Mark B3 fully resolved only once the response explicitly acknowledges the unresolved original cause and reports a verified supported comparison or gives the exclusion rationale requested by the reviewer.
5. Keep scType in the fair-comparison and paired-statistics work requested by R3 B2/B4. Original paper reporting is not a substitute for this required revision. GBM scType results are a different dataset and are unaffected by the PTC zero-T history.

## Sources inspected

- `paper/comments.md:2` — July 9, 2026 editor/reviewer source date.
- `paper/comments.md:63` — Reviewer 3 B3, SignacX zero-T challenge.
- `paper/comments.md:65` — Reviewer 3 B4, scType F1 0.9392 and accuracy 0.8823.
- `results/hvg_ptc_20260916_v1/ptc_paper_baseline/paper_Table2_literal_cells.json` — direct V16 DOCX table extraction.
- `results/hvg_ptc_20260916_v1/ptc_paper_baseline/original_writing_note/paper.md:113` — explicit method names and no-T statement.
- `results/hvg_ptc_20260916_v1/ptc_paper_baseline/paper_Table2_historical_source_comparison.csv` — historical metric reproduction.
- `results/hvg_ptc_20260916_v1/ptc_paper_baseline/paper_S3_literal_metrics.csv:58` — literal scType overall, group, sample counts.
- `results/hvg_ptc_20260916_v1/ptc_paper_baseline/paper_S3_literal_metrics.csv:72` — literal SignacX overall, group, sample zero counts.
- `results/hvg_ptc_20260916_v1/ptc_experiments/archived_comparators/format_normalized_name_ontology.csv:110` — scType native spellings, T calls, and separate NK calls.
- `results/hvg_ptc_20260916_v1/ptc_experiments/archived_comparators/format_normalized_name_ontology.csv:126` — populated original SignacX native output.
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_comparator_replay/ENDPOINT_CORRECTIONS.md` and `manifest.json` — completed cached replay and limits.
- `handoff/ptc_paper_baseline_20260916/reconcile_historical_metrics.py` — exact original endpoint reconstruction procedure.
- `handoff/paper_claim_validation_20260917/ptc_comparator_replay.py` — cached re-evaluation and label rules.
