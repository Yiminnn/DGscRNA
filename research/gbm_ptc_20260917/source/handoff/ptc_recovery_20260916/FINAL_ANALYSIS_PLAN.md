# Final PTC analysis specification

Written before final aggregate results are inspected. Extends the frozen controlled
protocol; it changes no fitting, annotation, label rule or condition availability.

The primary output is terminal DL/refinement at confidence 0.90. Marker-only calls
and threshold 0.70 are explicit ablations. The primary geometry feature comparison
is HVG2000 minus all, direct UMAP2/HDBSCAN15/15, full RNA scoring with
CellMarker_AllTissues/mean and constant group-selected RNA2000 DL features.
Clustering ARI and annotation concordance with S2 are historical concordance,
not accuracy against independent truth. Productive high-confidence TCR-positive
recall, detection yield, coverage and RNA support are complementary endpoints;
TCR absence is not an established negative cell identity.

Eight samples map to four patients. Average paired sample differences within each
patient before inference. Enumerate all 16 paired sign flips for two-sided tests
and all 256 size-four patient bootstrap resamples for descriptive 95% intervals.
The minimum nonzero two-sided exact p-value is 0.125. No cell-level significance
or seed-as-replicate inference is allowed. Report four patient values and sample
availability. Correct the five HVG-count comparisons within each path and endpoint
using Holm. Main endpoints and all additional method comparisons are transparent;
PTC is an exploratory reconstruction, not a confirmatory independent cohort.

Compare the six feature levels across all seven representations and three
clusterers. Separately quantify HVG-by-PCA preprocessing interactions, direct
UMAP/HDBSCAN versus each other default, and variation across five seeds. Flag the
one float64 GMM numerical recovery and give the affected condition with that
sample excluded; never silently drop it from the main census.

Batch comparisons use matched R single-sample controls and pooled NONE, CCA and
Harmony, with identical selected features, RNA scoring, DL input, UMAP and
clustering settings. CCA integrated-assay scoring/DL is a separate intervention.
Report both groups, paired patient effects, balanced shared-lineage mixing and
RNA-state/neighborhood preservation. Sample, tissue and patient are confounded;
mixing alone does not establish biological correctness.

Restore all 17 libraries and three cutoffs. Patient-held-out context selection
chooses a panel/cutoff using the other three patients' mean apparent productive
TCR binary F1, then reports the held-out patient's endpoints. Ties use library
name then cutoff order mean, none, p050. This is held-out panel selection on
transductively fitted predictions: expression embeddings and DL models were not
refitted without the held-out patient. It is not fully inductive validation.
Report selection versus fixed AllTissues, Thyroid, Pubmed and HPA contexts.

Use the corrected compact-name/encoding audit for archived competitors. Keep
literal S3 binary flags and native-name mappings separate. Do not combine old
competitors and new two-group DG-scRNA runs into a claimed fair method ranking.
Report 623 encoding-only SCINA differences, unresolved S3 TCR 5,435 and DG flag
3,901 cell discrepancies, original archived-output parity, missing historical
weights, exact raw counts/QC, and the dbscan index recovery evidence.

All statistics, scientific verification, figures and notebook execution use SLURM.
The final aggregate requires all 510 evaluation units, 11,696 terminal conditions,
49,776 sample-stage rows, complete feature/partition provenance, and independent
full-forward verification of every saved model cache used by these conditions.
