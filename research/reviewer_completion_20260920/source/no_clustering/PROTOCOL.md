# A2: GBM cell-wise seed replacement

Frozen before candidate fitting. This is an explicit replacement of cluster-DEG
seed construction, not removal of clustering with every other operation unchanged.
No PTC jobs are included. Scope is the original121 GBM samples, with HVG2000 and
HVG5000 DL inputs and all original RNA genes available for marker scoring.

For cell c and panel t, score = singleton_factor × sum(max(0, normalized RNA
expression[g,c])) over unique present marker genes, divided by the full original
panel length. singleton_factor=0.8 when panel length≤1, otherwise1. Missing genes
contribute zero; no cluster, DEG statistic or reference truth is used. Library is
CM2_glioma_other. The unique positive maximum is a seed candidate; a zero maximum
or tie is Undecided. Accept it when maximum ≥ λ × mean(maxima across all cells in
that sample). The primary λ is1. Sensitivities are0,0.5,1,1.5,2. Rejection produces
Undecided for the unchanged historical expression MLP to process.

Each arm keeps the original normalized selected-HVG DL input, 256/128 model,
training split/initialization seed42, ten epochs and original threshold rounding.
Terminal0.90 is primary;0.70 is a named secondary endpoint. Marker-only and final
outputs remain separate. No-op/all-Undecided/insufficient-training states remain
visible in metrics and figure captions. All cells remain in the denominator.

Selection uses the already frozen patient folds. Within each budget and cohort,
average sample terminal0.90 macro-F1 within training patients, then average equally
across training patients. Choose the λ with the highest training-patient mean;
ties choose nearest1, then smallerλ. All samples from each test patient use that
fold's choice. Test labels do not select λ. Primary97-sample and full121-sample
cohorts are evaluated separately. Fixedλ1 and selectedλ results are both retained.

All five candidates are reported, including low coverage and no-training cases.
Comparison to the existing four original routes uses the identical fixed marker,
mean cutoff and terminal endpoint, by patient and budget; no route is selected by
test outcomes. A2 differs in scoring mechanism as well as availability of cluster
DEGs. A favorable original-route contrast cannot establish that every possible
cell-wise/no-clustering method is inferior. The comparisons are retrospective.

Patient bootstrap intervals (2,000 replicates, seed42) condition on frozen fitted
predictions and selected configurations; they do not include re-fitting or
selection uncertainty. Two-sided paired Wilcoxon tests, if computed, use patient
macro-F1 differences. Holm correction spans all available predefined comparisons
within each cohort (2budgets ×4routes ×2selection modes=16); zero differences are
handled explicitly. Coverage effects are descriptive. No cell-level independence
is assumed.

The fitting process reads no author labels. Evaluation occurs only after immutable
initial and terminal predictions are present. Fitting and result cache parents
are pinned by source/protocol, native prep, cell, marker and DL input hashes.
Original-score/DL regression and a TKU4163 pilot gate the full SLURM array. The
array has at most one running task, four CPUs and32GB per task; failed/OOM jobs
remain incomplete until repaired. Original input/result files are read-only.
