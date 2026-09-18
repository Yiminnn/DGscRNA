# Bounded PTC follow-up after verified GBM delivery

This is an execution specification for the accepted NEXT_STEPS.md. It is written
before new PTC follow-up fitting. Earlier PTC outcomes are already known; this is
not a prospective registration. Existing results remain intact.

## Gate and scope

Every PTC scientific entry point requires SLURM and the verified
`GBM_full_summary/GBM_FULL_DELIVERED.json` gate. Code and scheduler preparation may
precede that gate. Reuse the archived all-eight-sample reconstruction and all 29
fresh/controlled PTC preparation units, each with 17 libraries, three density
cutoffs, four clustering routes and terminal DL. No independent Pu cohort or new
Darmanis fit is added. Historical missing weights do not block approximate
reproduction, but reconstructed weights are never described as recovered weights.

The two groups remain NMT = MT-1, MT-2, N-1, N-2 and TTU = TU-1, TU-2, T-1, T-2.
Patient1 = N-1/T-1; Patient2 = N-2/T-2; Patient3 = MT-1/TU-1;
Patient4 = MT-2/TU-2. The archived integrated fit has all eight samples; fresh
group fits are separate and must not be relabelled as the historical fit.

## Endpoints and selection

Use the frozen original label simplification and independently frozen name-only
strict-T rules. Recompute each patient's metrics from saved predictions, without
refitting the completed grid. Preserve the original paper broad-T/non-T F1 and
hard-call AUC against its original binary table. Separately report productive,
high-confidence TCR detection, T-positive F1, TCR-positive recall, detection yield
among predicted T cells, coverage, and Unknown-as-error accuracy/macro-F1.
Undetected TCR is not established non-T identity; detection-agreement statistics
cannot establish biological annotation accuracy. S2/S3 agreement is concordance.

Primary configuration selection maximizes the mean of the three training
patients' strict-T F1 against productive high-confidence TCR detection. All cells,
including Unknown, remain in the denominator. Fix this assay-concordance objective
before the new selection analysis; do not switch objectives to favour a route.
Leave out one patient at a time, excluding both samples from label-based selection.
Expression integration and pseudo-label refinement still use the full group:
these are transductive label-holdout results, not inductive unseen-patient fits.
Tie breaking uses the frozen configuration order. The 17-library roster and all
three cutoffs have equal opportunities within each compared route/feature budget.

Report fixed-original-marker versus training-patient-selected marker crossed with
initial versus terminal calls, using the same selected context for both stages.
Also retain stage-specific choices as a clearly separate sensitivity. Aggregate
patients equally, give all four paired differences, exact 16-sign-flip inference
and all 256 four-patient bootstrap resamples. Seeds and cells are not patients.
With four patients the smallest two-sided sign-flip p-value is 0.125; no significance
claim is manufactured from cell-level replication.

## Frozen original anchors and seed controls

The historical anchors are NMT CellMarker_Thyroid/none/PCA30-SNN and TTU
Pubmed_34663816/mean/UMAP2-HDBSCAN, on their archived all-eight-sample preparation.
Both are retained even if a different route ranks higher under a different endpoint.

New representation controls use fresh within-group CCA2000 and CCA5000, each with
PCA30-SNN and PCA30-to-UMAP2-HDBSCAN. This locks a competitive SNN comparison with
the same feature budget as UMAP-HDBSCAN; it is not chosen after new seed results.
Seeds are 0, 1, 2, 3, 42. Representation controls vary PCA and UMAP seeds and rerun
clustering, while keeping the CCA fit, feature identities, expression, SNN seed0,
density rules and MLP seed42 fixed. Existing seed42 predictions are the reference.
This tests conditional representation stability, not integration-seed stability.

MLP controls vary only model initialization across the same five seeds; known-cell
split seed42, architecture256/128, ten epochs, optimizer, cell order and batches
stay fixed. For each fresh configuration, retain the historical group marker plus
every distinct marker/cutoff selected by the four training-patient folds from the
completed seed42 grid. Freeze these choices before new seeds are fitted; the
held-out labels do not select a seed or change the chosen marker. Use the original
anchors as additional MLP controls. Native defaults must reproduce the completed
default predictions; numerical differences are documented before release.

## Uniform marker-retention mechanism control

Preserve the original geometry-only arms, including NMT's missing-CD3D limitation.
For each group use its completed all-gene CCA expression and a single scoring gene
set: fixed CCA2000 genes union all symbols in all 17 archived marker libraries,
intersected with that group's eligible all-gene CCA universe. The same scoring
genes and values, full original panel denominators and fixed 2000-gene DL matrix
are used for geometry2000, geometry5000 and geometryall, on PCA-SNN and UMAP-HDBSCAN.
This creates six conditions per group. It does not hand-insert a successful marker
or change labels. Missing markers outside the eligible universe remain missing.

Reuse and verify the exact completed geometry coordinates and cluster partitions;
recompute DEGs and density scores because the scoring assay has changed. All 17
libraries and three cutoffs reach terminal DL in these 12 conditions. Record
marker coverage and strict-T seed counts before and after the uniform rule. The
contrast is conditional on the all-gene CCA fit, not a replacement for CCA2000.

## Checks, resources and delivery

Default checks compare copied partitions, complete initial calls and terminal
labels/probabilities against the completed native-R reference. Perform the first
group preparation/scoring/terminal resource pilots before expanding. Use bounded
parallel SLURM jobs with immutable sources, actual resource accounting, checkpoint
hashes and explicit OOM/timeout recovery. Preserve failed attempts.

Every new clustering condition gets a figure on fixed, recorded display
coordinates. Extend the original notebook with the patient-holdout, module,
seed and marker-retention evidence; update the English decision tree and the
same authorized OneDrive directory. Retain original competing-method results with
their original input/reference conditions disclosed. The GBM matched-method
comparison is separate evidence, not proof that the PTC historical comparisons
were matched. Report original-winner, tuned-winner, ties and group-specific
counterexamples according to observed results.
