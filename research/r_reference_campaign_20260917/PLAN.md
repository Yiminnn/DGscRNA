# R-reference PTC ablation and reviewer datasets — 2026-09-17

The user explicitly requested continued PTC ablations and all previously selected
reviewer datasets with multiple biologically justified marker libraries. This
authorizes a new experimental phase; it does not resolve historical refitting
differences or the unreconciled manuscript Accuracy row. Preserve all old outputs.

## Scope and endpoints

- PTC: archived all-eight-sample checkpoint is the paper reference; NMT and TTU
  retain separate reporting and original selected routes. New within-group
  integration is an explicitly named intervention. All 17 original marker sets
  and all three original score cutoffs are retained.
- Reviewer roster: Baron human, Muraro, Segerstolpe, Xin, Immune_ALL_human, HCL,
  brain_GBM, breast_TNBC, colorectal, kidney_ccRCC, blood_DLBCL. This is the
  previously user-approved human cohort roster documented in
  `handoff/deck_datasets_provenance.md`, not every example in a cited review.
- No independent Pu cohort. GSE184362-derived markers remain included.
- Every annotation condition reaches terminal original-style DL/refinement, or
  reports its structural/no-op state. Marker-only calls are the no-DL ablation.
  Saved historic labels are concordance endpoints, not independent truth.

## Frozen comparisons

Use explicit R Seurat-v4-compatible DEG statistics, the archived density rule,
and the validated legacy PyTorch MLP implementation. No simplified-package
defaults are silently substituted. Keep cells and evaluation mapping fixed.
Record normalization, anchor genes, geometry genes, scoring genes, DL genes,
batch scope and random seeds separately for every condition.

1. Recover the complete original PTC marker × cutoff × four clustering grid on
   archived CCA geometry. Include the two reproduced selected routes as controls.
2. PTC ablation: HVG count, geometry, clustering, marker context, cutoff, DL and
   correction/scope. Isolated geometry contrasts hold scoring and DL inputs fixed;
   changing expression correction is a separate combined-workflow contrast.
   The 12 isolated-geometry controls use one all-shared-gene CCA fit per group,
   the same 2,000 scoring/DL genes and a byte-identical DL matrix; only the PCA
   feature budget varies (500/1k/2k/3k/5k/all). They are conditional on that fixed
   CCA fit, not replacements for the original CCA2000 reference. Together with
   17 combined preparation conditions and the archived baseline there are 30
   PTC analysis units; 70 reviewer analysis units give 100 total (including the
   explicitly paired colorectal CCA-feasibility sensitivity below).
3. Run reviewer datasets with the same reference algorithm and multiple marker
   contexts. Use available raw counts; explicitly describe already-normalized
   data when counts are unavailable. Do not pass RPKM as raw UMI counts silently.
4. Marker candidates are frozen from actual sampled organs and normal/disease
   context before evaluating predictions: primary tissue, relevant immune/stromal
   contexts, sampled metastatic/extranodal organs, unions and all-tissue controls.
   Preserve native normal/cancer and tissue prefixes; retain full denominators.
5. Hold evaluation labels out of annotation. Summarize the full marker grid, not
   a best-on-test claim. Where a selector is evaluated, select contexts using
   training donors and report held-out donor performance separately.

## Resources and delivery

All scientific data inspection, fitting, evaluation, tests and plotting run in
SLURM. Stage reusable expression, geometry and DEG outputs once; independent DL
conditions run in capped arrays. Pilot before expansion, inspect MaxRSS, preserve
failed logs, retry infrastructure errors without changing scientific conditions.
HCL requires a separately documented memory/scaling strategy after inventory.

Results: `results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/`.
Update the existing result notebook in place, retain all earlier GBM material,
show every clustering result and terminal annotation grid, and synchronize new
deliverables to the already authorized OneDrive directory. Keep statistical
claims limited to observed datasets/conditions; random seeds are not patients.

Primary CellMarker source: https://bio-bigdata.hrbmu.edu.cn/CellMarker2.0/index.html
and https://doi.org/10.1093/nar/gkac947. Freeze the cached human workbook checksum.

## Execution and interpretation details

- HCL retains all 599,926 cells across its 59 original tissue groups. Each group
  runs the same pipeline with its actual anatomical marker contexts; this is
  explicitly a tissue-conditional atlas analysis, not a pooled >100k clustering
  scalability experiment. Report its tissue count and aggregation rule.
- Fixed illustrative reviewer context is primary-normal CellMarker / mean;
  every other context remains in the complete grid. Descriptive maxima use the
  evaluation labels and must be named as such. Any donor-held-out label selector
  is transductive (all cells already entered unsupervised fitting), not an unseen
  donor refit. Preserve curated-semantic and common-lineage endpoints separately.
- Muraro uses the inspected raw.X integer matrix. Xin remains RPKM. The inspected
  Immune_ALL prepared X comes from its published mixed UMI/full-length count
  layer; fractional values are retained. Do not describe all inputs as raw UMI.
- All final model/NPZ/history hashes, exported prediction CSVs, known-label
  retention and 0.90/0.70 threshold reconstruction are audited. Preserve fixed
  PTC cell order, not merely the same cell set. Cross-condition feature identities
  and input hashes are tested explicitly in SLURM.
- Running R jobs use per-job immutable source copies after the recorded guard
  transition. The completed isolated-geometry preparation jobs are an exception:
  they ran before that driver received its guard. Their feature/cell-order/input
  audits establish output invariants, not exact historical execution bytes. The
  geometry transition record preserves only the source observed at the change.
  Earlier 16 reviewer units receive an independent unmodified
  original-density replay audit. This validates their outputs; it does not
  retroactively reconstruct transient source bytes. Later terminal jobs also
  freeze their driver and cap PyTorch at four threads, including high-memory
  all-gene tasks; the legacy MLP helper remains unchanged.
- Delivery includes a campaign-only SLURM task/step resource ledger and preserved
  failed-attempt logs. Active accounting rows remain explicitly non-final. A local
  delivery receipt alone is insufficient: its remote copy must be downloaded and
  hash-verified, and the finalizer must independently finish successfully in SLURM.

## Colorectal CCA-feasibility sensitivity (added after input audit)

The full 47,107-cell colorectal cohort includes donor HTA8_6004 with only three
curated cells. The full-cohort result remains an explicitly uncorrected RNA
fallback; it is not described as an executed CCA result. Add one CCA analysis
using a fixed donor-size threshold of at least 31 cells, needed for 30-dimensional
CCA: 28 donors and 47,104 cells. Preserve the excluded three cells in the original
cohort and record their IDs/reason. Reuse the exact 12 frozen colorectal marker
libraries, every clustering branch/cutoff and terminal DL. No labels or prediction
scores inform this exclusion. The two versions have different expression/scoring
universes and represent a combined workflow sensitivity, not an isolated batch
effect test. Count biological cohorts/cells only once in the roster.

## Quantitative correction diagnostics and completion checks

Six PTC correction controls (NMT/TTU × CCA2000, RNA, Harmony) receive exact k=30
neighbor diagnostics in saved PCA30/UMAP2. Queries are identical across arms,
sampled with seed42 within fixed archived population × sample strata. Batch
mixing is normalized by sample composition within each population and MT/N/TU/T
prefix; archived broad/native label purity, NMT tissue-state retention and TCR
detection homophily remain separate endpoints. Archived labels are concordance
references and sample/state confounding remains explicit. No fit is modified.

Final completion independently checks the frozen marker libraries × three cutoffs
× four named routes and all evaluation stages/scopes. New terminal task lists
are content-addressed and each arm has an exclusive writer lock; submitted job
registrations are persisted immediately. Stop/restart a dispatcher only while
holding its between-tick flock. Interrupted scoring caches are atomically
published or preserved under an invalid-cache name before recomputation.
