# Response to reviewers — working revision draft

Date: 2026-09-20. **Internal draft; not a submitted or completed response.**

This document covers all 35 requirements in the approved completion plan. It
preserves the original manuscript and the older response document. The response
paragraphs below are proposed wording. A statement about future revision is not
evidence that the final manuscript has been edited. Bracketed evidence slots must
be replaced with validated artifacts and final page/figure/table locations before
submission. Running jobs, completed jobs and scientifically accepted results are
distinct states. No new numerical result is claimed by this drafting task.

The reviewer source is `paper/comments.md`, dated 2026-07-09. Reviewer IDs follow
`handoff/reviewer_completion_plan_20260920/REVIEWER_COVERAGE.tsv`; the older July
reply uses a different subdivision and must not supply the current ID mapping.
The source manuscript is
`paper/submission_v16/DG_scRNA_04232026_V16_cell_report.docx`; source extraction,
hashes and supplement-container evidence are in `source_inventory.json`.

Fixed historical policy: preserve the original SignacX None / zero-predicted-T
result and original S3 labels. Preserve scType's actual result. Later nonzero
SignacX runs remain separate and are neither hidden nor forced to zero.

## R1-1 — Cost of the complete marker search

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:122–123`. Work packages: F,G.

**Proposed response**

We agree that the cost of one fitted workflow does not describe the cost of the complete selection procedure. The revision will separate shared preprocessing and integration, representation and clustering, differential expression, marker scoring, refinement, and aggregation. The accounting will identify unique computations, cache reuse, failed attempts, and retries, together with CPU-hours, elapsed time and peak memory. Existing single-workflow scalability measurements will remain a separate analysis.

**Pending acceptance slot**

[INSERT validated full-search resource table, hardware configuration, cache policy and wall-time definition; then bind exact manuscript and supplement locations.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `paper/comments.md:123`
- `handoff/paper_claim_validation_20260917/COMPLETION_EVIDENCE.json`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/terminal_execution_summary.csv`
- `results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/resources/`

## R1-2 — Interpretation of Unknown cells

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:124`. Work packages: D,G.

**Proposed response**

Unknown is an operational abstention category, rather than evidence of a new biological population. Existing GBM analyses examine expression and propagation of seed-label errors; the PTC extension will compare Unknown and annotated cells within supported broad classes and samples, including quality measures, marker modules and TCR support. Patient-level effects and uncertainty will accompany cell counts. Associations with doublet scores or intermediate expression will not be presented as proof of doublets or a transitional state.

**Pending acceptance slot**

[INSERT endpoint-matched PTC Unknown DE/QC results and gene-expression figures; retain any unsupported biological interpretation as a limitation.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/unknown_expression/all_patient_paired_gene_contrasts.csv.gz`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/unknown_summary/retained_seed_and_DL_new_errors.csv`
- `paper/comments.md:71`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/PTC_FOLLOWUP_REPORT_ZH.md`

## R1-3 — Framework figure numbering

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:125`. Work packages: G.

**Proposed response**

The framework illustration and all references to it will be reconciled in a preserved-copy revision of the manuscript. The main illustration will show the original preprocessing-to-terminal-annotation workflow; additional experimental branches will be placed in the supplement. Possession of a corrected figure is not equivalent to verification of the final submission PDF.

**Pending acceptance slot**

[INSERT final Figure 1 number, caption, manuscript pages and completed PDF cross-reference audit.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `handoff/paper_claim_validation_20260917/COMPLETION_REVIEW_20260919.md`
- `paper/comments.md:94`
- `paper/comments.md:113`

## R2-1 — External multiclass datasets

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:44`. Work packages: C,G.

**Proposed response**

We agree that a thyroid T-cell detection benchmark alone cannot establish broad annotation performance. Native-R outputs are already archived for the specified public-data roster, and the revision campaign will fill the remaining applicable comparator cells under declared input and information conditions. Multiclass results will include per-class and aggregate metrics, confusion matrices and Unknown coverage. Existing Darmanis results will be described with their label-compatibility limitations; they will not be used as a new tuning cohort.

**Pending acceptance slot**

[INSERT validated dataset-by-method coverage matrix and at least two complete external multiclass comparisons; cite all eleven predefined datasets without treating HCL tissue shards as one fitted job.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/summary/RESULTS_AND_INTERPRETATION.md`
- `results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/summary/cohort_sizes.csv`

## R2-2 — Contributions of clustering, marker selection and refinement

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:45`. Work packages: A,C,G.

**Proposed response**

The revised analysis separates the marker-selection and DL-refinement contributions and retains terminal annotation as the DG-scRNA endpoint. A no-clustering branch is being specified explicitly: cluster-DEG seed construction must be replaced by a cell-wise marker-seeding rule because the original density score is cluster-dependent. This contrast therefore estimates the effect of replacing seed construction, rather than removing clustering while leaving every other operation unchanged. Original-R anchors and paired representation comparisons will be reported separately.

**Pending acceptance slot**

[INSERT no-clustering protocol, anchor-parity evidence, matched terminal metrics, coverage and paired effects; do not insert results from the older simplified Python reply as native-R evidence.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/workflow_choice_summary/marker_DL_patient_paired.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/unknown_summary/retained_seed_and_DL_new_errors.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/marker_DL_paired_effects.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/controls_summary/representation_all_metrics.csv`
- `handoff/paper_claim_validation_20260917/prepare_representation_R.R`
- `paper/plan/IMPLEMENTATION_PLAN_CN.md:72`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/workflow_choice_summary/workflow_node_evidence.csv`

## R2-3 — Comparison with a GNN annotation method

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:46`. Work packages: C,G.

**Proposed response**

An existing scDeepSort comparison provides an actual GNN comparator; DG-scRNA itself is not redefined as a GNN to satisfy this request. The revision will document the pretrained reference, supported labels, missing classes and the target-data overlap audit for each usable run. Methods with different reference information will be reported in explicit information-condition groups, and unsupported malignant or other target classes will remain visible in the evaluation denominator.

**Pending acceptance slot**

[INSERT verified applicable GNN prediction rows, reference provenance and matched-class results; identify unsupported dataset-method cells rather than inventing predictions.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/comparison_summary/patient_heldout_summary.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/comparison_summary/diagnostics/`

## R2-4 — Expression profiles of unclassified cells

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:47`. Work packages: D,G.

**Proposed response**

The original code uses a terminal confidence threshold of 0.90; 0.70 is retained as a separately named sensitivity condition. The requested Unknown analysis will use the precise endpoint and threshold associated with each saved prediction. Differential expression and quality comparisons will be stratified by sample and supported broad cell class, with patient-level inference. Recovery strategies such as improving marker coverage will be discussed as proposed follow-up, not as already demonstrated recovery of novel cell types.

**Pending acceptance slot**

[INSERT endpoint-specific Unknown contrasts, QC/TCR associations and threshold labels; reconcile the manuscript's former 0.70-only statement.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/unknown_expression/all_patient_paired_gene_contrasts.csv.gz`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/unknown_summary/retained_seed_and_DL_new_errors.csv`
- `paper/comments.md:71`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/PTC_FOLLOWUP_REPORT_ZH.md`

## R2-5 — Sensitivity of parameters used by the implementation

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:48`. Work packages: A,B,G.

**Proposed response**

The original helper and terminal model define the parameters that can be meaningfully varied. Existing HDBSCAN, representation, MLP-width, epoch and seed controls will be reused. The new controls vary SNN neighborhood size and UMAP neighborhood size separately around their source defaults. The manuscript's alpha-weighted score and convergence-epsilon description are unsupported by the recovered implementation and will be corrected, rather than introducing new algorithm components solely to test them. Pilot findings will be identified as pilot findings.

**Pending acceptance slot**

[INSERT validated neighbor-control and learning-curve artifacts, actual parameter table, terminal outcomes and scope; no claim of a global optimum from a local sweep.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/GBM_full_summary/METHODS_CORRECTIONS.md`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/controls_summary/MLP_training_learning_curves.pdf`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/PTC_representation_seed_stability.pdf`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/PTC_MLP_seed_stability.pdf`

## R2-6 — Batch-associated annotation errors and biological structure

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:49`. Work packages: C,E,G.

**Proposed response**

The available CCA, no-correction and Harmony comparisons represent different complete workflows. They cannot establish a pure effect of moving an identical correction operation to another position. The revision will add comparable standard-tool predictions to the same PTC sample-level error and Unknown analysis, with mixing and biological preservation assessed within shared classes. Tissue and technical batch are partially confounded, so causal claims about a separately identified technical-batch effect will be avoided.

**Pending acceptance slot**

[INSERT same-endpoint DG/comparator batch-error tables for NMT and TTU, shared-class diagnostics and explicit confounding statement.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/summary/PTC_batch_biology_compact.csv`
- `results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/summary/RESULTS_AND_INTERPRETATION.md`
- `paper/plan/IMPLEMENTATION_PLAN_CN.md:14`

## R2-7 — Immunotherapy interpretation

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:50`. Work packages: D,G.

**Proposed response**

The biological discussion will connect supported T-cell expression states to relevant primary literature, while distinguishing expression association from functional activity, treatment response and validated predictive biomarkers. The cross-sectional cohort cannot by itself establish metastatic progression, a causal mechanism or therapeutic efficacy. Claims concerning therapeutic targets will therefore be framed as hypotheses requiring independent validation.

**Pending acceptance slot**

[INSERT verified primary references and the exact supporting marker/TCR panels; do not finalize a named target or predictive-biomarker claim without its source and cohort evidence.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `paper/comments.md:50`
- `paper/comments.md:71`
- `paper/comments.md:101`

## R2-8 — Distribution of marker-set performance

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:51`. Work packages: D,G.

**Proposed response**

The revision will show the distribution of terminal annotation performance across the predefined marker sources, along with their tissue, disease and assay provenance. Fixed-marker and training-patient-selected-marker analyses will be separated. Selected marker genes will be displayed in the same cells and terminal label endpoint used for interpretation. A library's name alone is insufficient to establish RNA-only provenance, and nominally mismatched libraries will not be omitted from the complete supplement.

**Pending acceptance slot**

[INSERT marker-source distribution, selected-gene dot/violin plots, gene-list table and selection record; reconcile actual candidate counts with the manuscript's historical 192 configurations.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/summary/marker_context_distributions.pdf`
- `results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/summary/marker_context_roster.csv`
- `paper/comments.md:89`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/comparison_summary/training_patient_choices.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/marker_evidence_summary/library_source_assay_coverage_audit.csv`

## R2-9 — Runtime and datasets larger than 100,000 cells

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:52`. Work packages: F,G.

**Proposed response**

The existing scalability campaign contains repeated independent CPU runs, including a largest single-run size of 120,000 cells. The revision will report this measured scope, hardware and memory accounting. It will not add separately fitted HCL tissue partitions and describe the sum as one large run. GPU memory is not applicable to the measured CPU implementation, and an unmeasured GPU advantage will not be claimed.

**Pending acceptance slot**

[INSERT verified 45-run resource table and curves, with supported methods, repeats and failure statuses; remove stale single-run 584,207-cell or GPU claims unless separately measured.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/scalability_summary/resource_repeated_mean_SD.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/scalability_summary/RESOURCE_INTERPRETATION.md`

## R3-A — Accurate positioning of graph and neural-network components

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:54–55`. Work packages: G.

**Proposed response**

The source audit identifies upstream clustering followed by context-dependent marker annotation and expression-MLP refinement. The terminal classifier receives an expression matrix rather than graph edges and does not implement message passing or joint graph representation learning. The revision will use this description consistently in the title framing, abstract, Methods and figures. The model is discriminative; the term generative annotation framework will be removed.

**Pending acceptance slot**

[Apply METHODS_REPLACEMENT.md to the preserved-copy manuscript and final figure captions, then verify the compiled document.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/GBM_full_summary/METHODS_CORRECTIONS.md`
- `paper/comments.md:54`

## R3-B1 — Multiclass evaluation versus T-cell detection

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:59`. Work packages: C,G.

**Proposed response**

TCR evidence supports a named binary T-cell detection endpoint; it is not complete multiclass ground truth or an independent validator of CD4, CD8 or Treg subtypes. The revision will separate this endpoint from public-data multiclass evaluation. Per-class, macro and weighted scores will use an explicit vocabulary, all eligible test cells and a declared Unknown rule. Historical Table 2 F1 and hard-call AUC will retain their recovered definitions rather than being relabeled as modern multiclass metrics.

**Pending acceptance slot**

[INSERT frozen evaluation vocabularies and multiclass tables/confusion matrices; keep the unresolved historical Accuracy definition explicitly labeled.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_comparator_replay/ENDPOINT_CORRECTIONS.md`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/original_paper_endpoint_reconciliation.csv`
- `results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/summary/RESULTS_AND_INTERPRETATION.md`
- `results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/summary/cohort_sizes.csv`

## R3-B2 — Leakage control and fair marker selection

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:61`. Work packages: C,G.

**Proposed response**

Selecting a marker set using all evaluation labels can inflate performance. The revision protocol keeps all samples from a patient in one fold and uses training patients for marker and parameter selection, with test patients used only for evaluation. Fixed-marker and selected-marker comparisons will be shown separately, and marker versus pretrained-reference methods will disclose their information conditions. Existing label-held-out but jointly integrated analyses are transductive and will not be described as fully inductive unseen-patient projection.

**Pending acceptance slot**

[INSERT fold manifests, selection logs, allowed marker/parameter budgets and compatibility matrix; apply descriptive-only labels where donor/study independence is unavailable.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/comparison_summary/training_patient_choices.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/heldout_workflow_patient_metrics.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/comparison_summary/patient_heldout_summary.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/comparison_summary/equal24_DG_scType_comparison.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_comparator_replay/ENDPOINT_CORRECTIONS.md`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_comparator_replay/paired_patient_gains.csv`

## R3-B3 — Historical SignacX zero-T result

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:63`. Work packages: C,G.

**Proposed response**

The original paper's SignacX entries are preserved as reported: None, with zero predicted T cells under the archived endpoint. This is not scType, whose original results remain unchanged. The recovered SignacX vector contains non-T predictions, so zero T calls do not imply a crash or absence of output. The original pooled integrated-input SignacFast run differs from the later per-sample raw-input Signac run in several ways; the causal reason for the historical zero is not isolated. Later cached nonzero predictions and any new comparable evaluation will be reported separately and will not overwrite, or be forced to reproduce, the historical zero.

**Pending acceptance slot**

[INSERT original input/native-label provenance appendix and a supported fair SignacX comparison or explicit exclusion rationale; do not use the historical zero alone as evidence of general superiority.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_comparator_replay/ENDPOINT_CORRECTIONS.md`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_comparator_replay/paired_patient_gains.csv`

## R3-B4 — Effect sizes and uncertainty against every comparator

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:65,116`. Work packages: C,G.

**Proposed response**

The revision will report absolute and relative differences against every eligible comparator rather than highlighting the largest relative gain over the weakest baseline. Paired analysis will follow the patient unit for GBM and display all four PTC patients, with appropriately limited exact paired inference and multiplicity control where used. Historical binary-call AUC will not be treated as a continuous-probability AUC. Existing statistical outputs are reusable only when their endpoint, predictions and information conditions match the declared comparison.

**Pending acceptance slot**

[INSERT complete matched-comparator effect and interval tables, tests, Holm family definition and support counts for new arms.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/comparison_summary/paired_patient_comparisons.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/paired_patient_contrasts.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/comparison_summary/equal24_DG_scType_comparison.csv`

## R3-C1 — Repeated runtime and memory measurements

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:69`. Work packages: F,G.

**Proposed response**

Runtime and memory reporting will distinguish independently repeated single-workflow scaling from the full-search ledger. The existing CPU scaling outputs will be presented with hardware, allocated resources, wall-time definitions, cache state and repeat variability. Maximum resident memory must be interpreted in relation to the measured process or SLURM step and cannot automatically be treated as the summed memory of every child process. GPU measurements are not available for the CPU-only campaign.

**Pending acceptance slot**

[INSERT resource definitions and verified repeated-run figure; link full-search accounting from R1-1 without conflating the two experiments.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/scalability_summary/resource_repeated_mean_SD.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/scalability_summary/RESOURCE_INTERPRETATION.md`

## R3-C2 — Evidence for subtypes and patient-level composition

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:71`. Work packages: D,G.

**Proposed response**

The revised interpretation will preserve the distinct original S2 and S3 annotation endpoints and display marker evidence and composition separately for each. Sample-level views and paired patient summaries will accompany pooled composition. The N-to-T and MT-to-TU contrasts each involve only two patients, not eight independent biological replicates. Marker expression and TCR support provide supporting evidence, but not independent subtype ground truth or proof of causal progression or a newly discovered population.

**Pending acceptance slot**

[INSERT endpoint-matched selected-marker panels, all-patient composition and paired displays, and revise unsupported causal or novelty language.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `paper/comments.md:71`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/PTC_FOLLOWUP_REPORT_ZH.md`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/original_anchor_patient_metrics.csv`
- `paper/comments.md:50`
- `paper/comments.md:101`

## R3-C3 — Refinement gains and propagation of seed-label errors

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:73`. Work packages: C,D,G.

**Proposed response**

The existing marker-by-refinement contrasts and error-transition analyses distinguish retained seed labels, newly assigned cells, corrections and new errors. The original refinement preserves initially assigned labels and fills only Undecided cells, so it cannot automatically correct every wrong seed. The revision will report benefits and harms conditional on the marker context and usable training state. No-op and invalid branches will be distinguished from fitted models.

**Pending acceptance slot**

[INSERT saved transition tables and matched coverage/performance figures for each retained endpoint; integrate new valid branches without dropping unfavorable transitions.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/workflow_choice_summary/marker_DL_patient_paired.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/unknown_summary/retained_seed_and_DL_new_errors.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/marker_DL_paired_effects.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/unknown_expression/all_patient_paired_gene_contrasts.csv.gz`

## R3-S1 — Clustering metrics and genuine cross-tool agreement

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:84`. Work packages: F,G.

**Proposed response**

Partition metrics and terminal annotation metrics answer different questions and will be labeled separately. The added agreement analysis will compare predictions from different tools on the same cells and vocabulary, including pairwise agreement, FMI and all-tool agreement where applicable. Agreement and disagreement strata will be evaluated against the declared reference, with Unknown treated explicitly. Agreement among tools is not evidence of correctness by itself, and voting among reference sets within one method is not a cross-tool consensus.

**Pending acceptance slot**

[INSERT validated cross-tool matrices, truth-stratified agreement counts and denominators; preserve separate partition ARI/NMI/FMI labels.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `paper/comments.md:84`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/comparison_summary/`
- `results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/summary/all_clustering_metrics.csv`

## R3-S2 — Parameters and unsupported formula terms

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:85,115`. Work packages: B,G.

**Proposed response**

The source-supported parameter table uses 2,000 variable or integration features, PCA30, the stated SNN and UMAP neighborhood sizes, R HDBSCAN minPts50, a 256/128 expression MLP and a primary terminal threshold of 0.90. The number 3,000 in the saved ScaleData command is min.cells.to.block, not the HVG budget. The source implements a DEG log-fold-change density score with a singleton penalty, not the manuscript's alpha mixture or epsilon iteration. These discrepancies will be corrected directly in Methods.

**Pending acceptance slot**

[Insert METHODS_REPLACEMENT.md and the actual sensitivity table; identify 3,000-HVG and 0.70-threshold conditions as named comparisons.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/GBM_full_summary/METHODS_CORRECTIONS.md`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/controls_summary/MLP_training_learning_curves.pdf`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/PTC_representation_seed_stability.pdf`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/PTC_MLP_seed_stability.pdf`

## R3-S3 — Normalization at each stage

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:86`. Work packages: G.

**Proposed response**

The recovered reference helper uses Seurat LogNormalize with a scale factor of 10,000. It does not establish a TMM/edgeR normalization stage for the reproduced annotation workflow. The revised text will describe normalization, feature selection, sample-level scaling, integration and downstream expression inputs separately, including the distinction between the RNA and integrated assays. A different normalization procedure will not be retrofitted to make the code match the earlier prose.

**Pending acceptance slot**

[Replace V16 paragraphs P0600 and P0664 and verify all Methods/figure references; report any separate analysis's normalization under its own endpoint.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/GBM_full_summary/METHODS_CORRECTIONS.md`

## R3-S4 — Training duration and validation curves

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:87`. Work packages: B,G.

**Proposed response**

The original ten-epoch endpoint will remain the reproduction anchor. The diagnostic campaign records training and held-out pseudo-label validation loss and accuracy at every epoch, with prespecified 5/10/20/30 checkpoints and seeds. This internal validation tests generalization to marker-derived labels and is not independent biological truth. Test labels will not be used for early stopping, and a plateau or optimum will not be asserted before the actual curves are inspected.

**Pending acceptance slot**

[INSERT validated learning curves, original-endpoint parity and checkpoint outcomes; record legitimate no-training cases instead of fabricated curves.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/controls_summary/MLP_training_learning_curves.pdf`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/PTC_representation_seed_stability.pdf`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/PTC_MLP_seed_stability.pdf`

## R3-S5 — Incomplete TCR capture

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:88`. Work packages: D,G.

**Proposed response**

Absence of a productive TCR observation does not prove a cell is non-T. The revision will present the existing TCR-definition and T-core threshold sensitivities, including patient-level proportions of transcriptionally supported T cells lacking a productive receptor call. This is a definition-dependent proxy for missing TCR evidence, not an identifiable true capture false-negative rate without independent T-cell truth. The distinction will be stated alongside the binary reference assumption.

**Pending acceptance slot**

[INSERT definition table, numerator/denominator and per-patient sensitivity figures; retain limitations for negative-reference labels and subtype validation.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `paper/comments.md:88`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/original_anchor_patient_metrics.csv`

## R3-S6 — Expression of the selected context markers

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:89`. Work packages: D,G.

**Proposed response**

The revised marker figures will show actual gene expression, rather than only marker-score distributions. Each panel will identify the marker source, gene list, selection rule, assay, saved run and terminal labels. The historical NMT Thyroid/PCA-SNN/none and TTU Pubmed/UMAP-HDBSCAN/mean contexts will remain distinct. Missing or unexpressed markers will be visible rather than silently removed from the provenance record.

**Pending acceptance slot**

[INSERT validated selected-context dotplots/violins and full context-to-gene tables; GBM-first execution does not mean PTC panels are already complete.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `paper/comments.md:89`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/comparison_summary/training_patient_choices.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/marker_evidence_summary/library_source_assay_coverage_audit.csv`

## R3-S7 — Per-class and macro-averaged metrics

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:90`. Work packages: C,G.

**Proposed response**

The revision will provide per-class precision, recall and F1, together with macro and weighted summaries, confusion matrices and coverage. The class universe and handling of absent, unsupported and Unknown labels will be explicit. All eligible test cells will remain in the denominator under the prespecified scoring rule. The historical PTC non-T-positive binary F1 is retained only in the historical reconstruction table and is not renamed macro-F1 or T-positive F1.

**Pending acceptance slot**

[INSERT endpoint-matched tables for all new valid method/dataset rows and verify averaging definitions in captions.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `handoff/paper_claim_validation_20260917/DARMANIS_SCOPE_UPDATE_ZH.md`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/protocol/`
- `handoff/meeting_0804/ACTION_ITEMS_0804.md`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/summary/primary_fixed_marker_24.csv`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/summary/aggregate_manifest.json`

## R3-M1 — Figure references and classifier dimensions

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:94–95`. Work packages: G.

**Proposed response**

The source and main model description specify hidden-layer widths of 256 and 128. The revision will reconcile this with the conflicting network schematic and correct framework figure references throughout the document. The main figure will include preprocessing and the terminal classifier with an expression input; it will not depict nonexistent graph-message-passing layers.

**Pending acceptance slot**

[INSERT final consistent figure/caption and audit the compiled manuscript for unique numbering and resolvable references.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `handoff/paper_claim_validation_20260917/COMPLETION_REVIEW_20260919.md`
- `paper/comments.md:94`
- `paper/comments.md:113`
- `paper/comments.md:95`
- `paper/comments.md:99`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/GBM_full_summary/METHODS_CORRECTIONS.md`

## R3-M2 — Cell attrition from delivery to analysis

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:96`. Work packages: G,H.

**Proposed response**

The original delivered GEX files contain 110,497 barcodes, whereas the historical downstream analysis endpoint contains 92,404 cells. These are different stages. The manuscript's statement that 110,497 cells were verified as retained after the complete analysis QC is not supported by the delivery record. The revision will provide a sample/barcode membership ledger for each available checkpoint. Missing intermediate reasons will be labeled unavailable rather than assigned to QC or doublet filtering by subtraction.

**Pending acceptance slot**

[INSERT cell-level stage ledger and validated per-sample totals; use SOURCE_GAPS.md for the known delivery counts and unresolved filtering provenance.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/geo_submission_v1/PTC_GEO_submission_20260916/README.md`
- `paper/comments.md:96`

## R3-M3 — Wet-laboratory and historical software provenance

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:97–98`. Work packages: H,G.

**Proposed response**

The source documents contain incompatible viability thresholds and strainer pore sizes. We cannot select a value on the basis of the current computational environment. Six original GEX reports identify Single Cell 5-prime PE chemistry, conflicting with the manuscript's 3-prime wording; two GEX reports are unavailable. The historical DoubletFinder version also remains unverified: a version installed for later recovery is not proof of the original run. These values require the original laboratory or run records before the final Methods can use one definitive specification.

**Pending acceptance slot**

[Resolve SOURCE_GAPS.md with original viability/protocol/kit records and original session or environment evidence; additionally reconcile the reported Cell Ranger versions.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `paper/comments.md:97`
- `paper/comments.md:98`
- `results/hvg_ptc_20260916_v1/geo_submission_v1/PTC_GEO_submission_20260916/README.md`

## R3-M4 — Software citations and Key Resources Table

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:98–99`. Work packages: G.

**Proposed response**

The revised resources table will align each tool, version, URL and primary citation. The original table attributes DoubletFinder to Xun rather than McGinnis and links the SignacX row to the Signac chromatin-analysis citation; these are distinct sources. The normalization and clustering package entries also need to match the actual original-R implementation. The final audit will check SCINA, scType, scCATCH/CellMatch and all remaining citations against primary sources.

**Pending acceptance slot**

[INSERT complete source-verified bibliography/resource mapping and final text cross-reference check; no new external citation verification is claimed by this local drafting task.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `paper/comments.md:95`
- `paper/comments.md:99`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/GBM_full_summary/METHODS_CORRECTIONS.md`

## R3-M5 — Claims supported by the evidence

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:100–102`. Work packages: G.

**Proposed response**

The revision will describe a discriminative expression classifier following clustering and context-dependent marker annotation. Claims of universal superiority, a generative model, independent discovery of novel populations or causal metastatic progression will be limited to what the data establish. New comparisons will be reported for their explicit datasets, candidate workflows and endpoints; original-R provenance alone does not prove that every step is optimal. All valid comparisons will remain available in the supplement.

**Pending acceptance slot**

[Apply a sentence-level claim audit after the validated results arrive; bind each retained performance or biological claim to its figure/table and scope.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/GBM_full_summary/METHODS_CORRECTIONS.md`
- `paper/comments.md:54`
- `paper/comments.md:71`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/PTC_FOLLOWUP_REPORT_ZH.md`
- `results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/summary/RESULTS_AND_INTERPRETATION.md`
- `results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/summary/cohort_sizes.csv`

## R3-T1 — GEO accession and reviewer access

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:112`. Work packages: H.

**Proposed response**

A local GEO preparation package is available with processed files, sample crosswalks and provenance. It is not an assigned GEO accession or verified reviewer-access service. The currently held package records missing raw sequencing files and unresolved metadata. These gaps will be resolved through the original data providers and an authorized submission before a data-availability statement claims public or reviewer access.

**Pending acceptance slot**

[INSERT actual accession, checked reviewer access, raw/processed file status and release information only after authorized submission; local package creation does not close this item.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `results/hvg_ptc_20260916_v1/geo_submission_v1/PTC_GEO_submission_20260916/README.md`
- `paper/comments.md:112`

## R3-T2 — Accessible original supplementary material

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:113`. Work packages: G,H.

**Proposed response**

The original Figure S1 and Tables S1-S4 are locally present and are preserved. The four XLSX containers pass CRC checks and expose readable workbook metadata, but this is not a full endpoint/content or final reader test. The revised supplement index will retain the original files, provide legible figure/table previews and relate them to the new Supplement A/B analyses. Historical S2 and S3 will remain separate annotation endpoints.

**Pending acceptance slot**

[Complete SUPPLEMENT_INVENTORY.md checks, correct the S1 workbook's misleading internal sheet title in a revision copy if appropriate, and validate all final links/captions and downloadable tables.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `handoff/paper_claim_validation_20260917/COMPLETION_REVIEW_20260919.md`
- `paper/comments.md:94`
- `paper/comments.md:113`

## R3-T3 — Reproducible examples and permanent archive

Reviewer source: `/fs/scratch/PCON0080/yimin/dgscrna/paper/comments.md:114`. Work packages: H,G.

**Proposed response**

The reproducibility package will distinguish replay of saved original Sup endpoints from an original-R refit whose numerical tolerance is explicitly defined. The former cannot restore unavailable model weights, and neither can manufacture the unresolved historical Accuracy formula. The package will include frozen inputs and environments, portable entry points, expected outputs and SLURM instructions. A clean-environment execution and an actual permanent archive identifier remain separate acceptance requirements.

**Pending acceptance slot**

[INSERT successful fresh-environment Table 2/GBM example logs and immutable release metadata; add DOI only after authorized publication and access verification.]

**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.

- `handoff/paper_claim_validation_20260917/README.md`
- `handoff/paper_claim_validation_20260917/COMPLETION_REVIEW_20260919.md`
- `paper/comments.md:114`
- `handoff/paper_claim_validation_20260917/COMPLETION_EVIDENCE.json`
- `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/PTC_summary/original_paper_endpoint_reconciliation.csv`

## Finalization gate

Before any response is described as complete, replace every acceptance slot with
verified evidence or an explicit source-supported limitation and the corresponding
manuscript change. Check all35 IDs against the approved acceptance table, including
external-source requirements for wet-lab facts, raw reads, actual GEO access and
permanent archival. Produce and inspect the final submission PDF and supplements.
Do not convert partial computational coverage into a claim that all reviewer
requirements have been closed.
