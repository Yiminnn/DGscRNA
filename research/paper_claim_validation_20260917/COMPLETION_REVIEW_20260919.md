# Completion review: accepted GBM and PTC follow-up plan

The computational and delivery scope in `NEXT_STEPS.md`, with the later
`DARMANIS_SCOPE_UPDATE_ZH.md` and accepted approximate PTC reproduction standard,
is complete. This review reads the actual saved results, producer-specific
completion flags, SLURM accounting, archive audits and remote receipts. It does
not rerun an experiment or repeat a completed package audit.

Paths below refer to the existing workspace/OneDrive result tree. The result root
is `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917` (abbreviated `O`).
The code archive does not contain raw data, model weights or rendered notebooks.
Exact receipt links and current-file hashes are in `COMPLETION_EVIDENCE.json`.

| Accepted requirement | Completed scope and authoritative evidence |
|---|---|
| Native R, SLURM execution and final DL/refinement | Frozen execution sources, job-linked manifests and native-R parity records are retained. Marker-only results are explicitly ablations. Terminal execution states distinguish training, cache reuse, no-op and untrainable cases. |
| GBM cohorts, pilots and complete six-budget/four-route grid | Fixed 121-sample/59-patient cohort and 97-sample/55-patient primary analysis; 726 units and 2,904 clustering results. `O/verification/pilot_manifest.json`, `O/summary/aggregate_manifest.json`, `AGGREGATE_COMPLETE`, and the completed core archive audit cover inputs, terminal results and figures. |
| Feature, geometry, marker, DL and parameter contributions | 363 geometry controls, 54 MLP controls and 180 representation controls, with learning curves and the exact 2,000-feature anchor. See `O/controls_summary/manifest.json`, `O/workflow_choice_summary/manifest.json` and their tables. The latter contains 21 workflow-node evidence records. Limited parameter pilots are not represented as a global parameter search. |
| Marker sources and patient-held-out comparisons | All 16 GBM libraries are retained; 13 eligible libraries enter primary selection. Source/label overlap is disclosed. Five comparator methods have 121 evaluated sample outputs each; fixed-partition and equal-24-configuration comparisons are separated. See `O/marker_evidence_summary`, `O/comparison_summary` and the GBM follow-up archive audit. |
| Unknown and error propagation | Both Unknown diagnostics and patient-paired expression analysis are complete for the 97 primary samples; 26,590 supported gene contrasts are reported with the original support criteria. These associations are not claimed to establish new cell types, doublets or therapeutic targets. |
| Resource experiments | All 45 runs are complete: three methods, five sizes (10k, 30k, 50k, 100k and 120k), three independent processes per condition. Actual timing, job/step memory and allocation details are retained in `O/scalability_summary`. OS caches are not claimed controlled. |
| Finish GBM before new PTC computation | `O/GBM_full_summary/GBM_FULL_DELIVERED.json` opened the PTC gate at 2026-09-19 11:28:44 UTC, after the full GBM download verification. PTC execution began after that gate. Core and follow-up GBM delivery receipts and completed content audits remain intact. |
| PTC baseline and fixed groups | Original all-eight-sample CCA and the requested NMT/TTU group integrations remain separate. The original NMT SNN/CellMarker_Thyroid and TTU UMAP-HDBSCAN/Pubmed_34663816 anchors are retained. Thirty original grid units were re-evaluated; native scoring and reconstructed terminal parity passed. Paper non-T F1/AUC, strict-T/productive-TCR concordance and S2/S3 agreement have distinct definitions. |
| PTC selection, five seeds, four-way analysis and marker retention | All 22 control units, 44 clustering routes and 50 MLP task groups completed. The original 17 libraries and three cutoffs are retained. All 12 marker-retention route audits passed fixed-partition, geometry, DL-input and shared-gene DEG checks. `O/PTC_summary` records four-patient label holdout, marker × DL effects, initialization versus representation variation, and paired contrasts. |
| PTC terminal accounting | All 811 requested follow-up conditions are accounted for once: 667 fresh training executions, 123 exact cached-terminal reuses, 19 uncached all-known no-ops and two uncached no-known-label conditions. There are no invalid terminal results. Raw single-known-class and cached training states remain separately visible in the ledger. The two archived seed42 parity pilots count once within the 50 MLP task groups. |
| Cached PTC competitors and endpoint correction | `O/PTC_comparator_replay` joins saved predictions to all 92,404 evaluation cells, preserves patient differences and Unknown-aware denominators, and separates SignacX explicit-T CellStates from coarse TNK CellTypes. No SignacX model was refitted; the original no-lymphoid output is not mislabeled as a crash. |
| Existing notebook, all clustering figures and English workflow | The canonical file remains `notebooks/dgscrna_results.ipynb`: 801 cells, including all 694 prior cells unchanged and 107 new cells, with zero execution errors. Its pre-PTC backup is delivered. The 44 new PNG/PDF pairs, existing 2,904 GBM clustering results, and English decision tree remain indexed in the original directory. Five new PTC summary plots and representative NMT/TTU cluster plots were visually inspected. |
| Same OneDrive directory, no replacement site | The existing `onedrive:work_od/share/dgscrna_GSE274546_TKU3186/01_report/GBM_PTC_results_20260916` contains the organized results and updated index. All 31 local targets in its main README resolve. The PTC main package has 388 download-verified files; its content audit verified 105 archives, all 811 requested terminal conditions and both notebook versions. The audit supplement was independently copied, download-checked and its receipt read back. |
| Reviewer cohorts and exclusions | Existing 11-dataset native-R results remain available through the previous campaign index. HCL's 599,926 cells across 59 tissue units are not called a single 600k fit. No new Darmanis tuning or independent Pu cohort was added; the original GSE184362-derived marker library remains in scope. Old/raw results were preserved. |
| Repository and account scope | Scientific execution code was already pushed on `align-r-reference`; final closure records are archived only in `research/paper_claim_validation_20260917`. The verified local GitHub account is `Yiminnn`. Other branches, the Python package and the unrelated untracked `requirements-lock.txt` are outside this change. No unknown connector or Remote Control was used. |

The final PTC job `7398160` completed with exit `0:0` at 2026-09-19 14:59:56 UTC.
Its main download verification, first package-content audit and audit-supplement
verification all completed in that allocation. `O/PTC_summary/PTC_FULL_DELIVERED.json`
links those actual receipts; it is not a job-submission or fitting-only flag.

## Scientific conclusion retained in the delivered reports

Completion does not establish the desired universal optimality claim. In GBM,
equally selected workflows yielded primary patient-mean macro-F1 0.284639 for
PCA/SNN and 0.232700 for UMAP/HDBSCAN. This compares selected workflows, including
their feature/marker/cutoff choices, rather than an isolated clustering effect.
The geometry-only and DL-only controls, coverage changes, competitor results and
counterexamples remain reported.

PTC retains its different original branches. Four-patient mean strict-T versus
productive-TCR F1 was 0.906657 for the reconstructed NMT anchor versus 0.871218
for the training-patient-selected fresh workflow, and 0.587918 versus 0.820465
for TTU. These are assay-concordance endpoints, not multiclass annotation accuracy.
TCR absence does not prove non-T identity. Uniform marker retention restores the
previous zero NMT scores while keeping geometry and DL inputs fixed, so the
earlier zeros cannot establish an inferior geometry choice. Original and retained
scoring universes are both preserved.

Only four paired patients are available in PTC; the minimum exact two-sided
sign-flip p-value is 0.125. Seeds and cells do not increase biological replication.
Selection is retrospective and integration/refinement remain transductive.
Historical trained weights and the source of the original Accuracy row remain
unavailable, as already accepted and disclosed. Wet-lab validation, permanent
repository DOIs and manuscript submission are not claimed completed by this
computational scope. No required experiment or results-delivery task remains.
