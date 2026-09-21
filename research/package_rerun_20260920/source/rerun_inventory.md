# Packaged R-reference rerun inventory

Inventory date: 2026-09-20 EDT. This is a read-only metadata inventory and rerun plan,
not a job submission or package acceptance. No expression/prediction data were
loaded, no jobs stopped, and no remote service accessed.

## Immediate launch sequence

1. Freeze the built package wheel/source commit and dependency lockfiles. Install
   the built artifact into a fresh prefix; invoking source directly is not enough
   to prove that the published package contains its R/DL resources.
2. Fresh packaged TKU4163/HVG2000/UMAP2_HDBSCAN_R/CM2_glioma_other/mean run.
   Compare cells/features, normalized DL matrix, PCA/UMAP, partitions, initial
   calls, split, complete probabilities and terminal090/070 with the existing
   native-R reference. Reuse zero old models or predictions. RDS serialization
   bytes can differ after relocation; compare scientifically relevant fields.
3. Fresh NL022 (3,628 cells) as the non-fixture input; cover HVG2000 and all genes.
   Also run SN040 (10,135 cells) as the large-memory pilot. TKU4163 has 178 cells.
   At least one pilot must cover all four routes and every marker/cutoff arm,
   including explicit no-training terminal states. A single L00_mean success
   cannot certify the complete 48-arm scorer.
4. Only after the above parity and resource checks, release 726 preparation jobs
   from `rerun_inventory.core_tasks.json` through the **actual new package API**.
   Its rows specify inputs/configurations, not speculative CLI flags. Output in
   `results/hvg_ptc_20260916_v1/package_reference_rerun_20260920/` or a parent-chosen
   fresh version; never reuse the old result root or old DL cache.
5. Independently verify all terminal states, re-evaluate using the frozen Lfine
   endpoint, and refresh the existing compact notebook locally. Then progress
   through the GBM controls/extensions, PTC baseline gate, PTC/public experiments.

The current input counts already contain author-retained GBM cells. Preserve this
roster and the >=3-cell gene filter for a comparable rerun. Applying a new raw-QC
or DoubletFinder recipe would be a new experimental condition, not parity.

## Experiment families

| Priority | Family | Frozen scope | Existing state / rerun treatment |
|---|---|---|---|
| 0 | Installed-artifact parity | TKU4163, NL022, SN040; HVG2000 plus all; full routes/arms on at least one pilot | Fresh fits required; fixed TKU4163 example already had clean-environment parity, but a new package needs its own proof. |
| 1 | GBM core | 121 samples × HVG500/1000/2000/3000/5000/all × 4 routes × 16 libraries × none/mean/0.5 | 726 preparations, 2,904 partitions, 139,392 terminal configurations already complete natively. Rerun in the package with isolated fresh caches. |
| 2 | Geometry-only controls | 121 × HVG2000/5000/all × 4 routes × 2 fixed markers/mean; fixed RNA scoring and normalized HVG2000 DL | 363 preparation/control tasks; 2,904 terminal configurations. Reuse newly packaged core geometry, hold DL/scoring fixed exactly. |
| 2 | MLP-only controls | 3 pilots × 2 budgets × 9 model settings × 4 routes × 2 markers | 54 task groups, 432 terminal configurations. Seeds0/1/2/3/42, widths128/64 and512/256 versus256/128, epochs5/10/20; split42 fixed. |
| 2 | Representation parameters | 3 pilots × 2 budgets × 30 geometry/parameter settings × 2 markers | 180 tasks/360 terminal configurations; UMAP10/30, noDR, minPts25/100, resolution0.25/1 and representation seeds0/1/2/3. |
| 2 | Reviewer B neighbors/learning | 66 neighbor conditions × 2 markers; 18 learning trajectories with epochs5/10/20/30 | Native fits complete; 408 terminal threshold rows independently Lfine-verified. Port actual controls, not just saved metrics. |
| 2 | Reviewer A2 no-cluster seeds | 121 × HVG2000/5000 × λ0/0.5/1/1.5/2 | 242 sample/budget units, 1,210 terminal conditions. This replaces seed construction; it is not an unqualified deletion of clustering. Fits complete; Lfine replay has 6,534 stage rows including original-route anchors. |
| 2 | Reviewer A1 seven-space comparison | 121 × HVG2000/5000 × noDR/PCA2/FA2/ICA2/Isomap2/UMAP2/TSNE2 × KMeans6K/GMM6K/HDBSCAN | 1,694 representation units, 22,022 candidate partitions, 66,066 initial/090/070 metric rows. Running now in another research controller; preserve it and stage any package rerun separately. |
| 3 | GBM method comparators | DG, scType, scCATCH, SCINA, SingleR, corrected scDeepSort | Existing native/corrected outputs available; new packaged DG should be re-compared against fixed comparator predictions where their input contract is unchanged. Do not retrain unrelated methods merely to relabel unchanged results. |
| 4 | PTC historical baseline and grouped baseline | Fixed S2/S3 endpoints; NMT Thyroid/none/PCA30-SNN; TTU Pubmed_34663816/mean/UMAP2-HDBSCAN | Reconcile after packaged GBM acceptance. Historical all-eight CCA is distinct from fresh grouped CCA. Never promise lost historical weights will reproduce cell-exact predictions. |
| 5 | PTC existing grid | Archived baseline1 + combined preparations17 + isolated geometry12; each4routes×17libraries×3cutoffs | 30 units/6,120 terminal conditions complete. Fresh two-group CCA budgets500/1k/2k/3k/5k/all, NONE/Harmony controls; no independent Pu. |
| 5 | PTC followups | Representation/retention22 preparations,44 route scores,696 terminal arms;50 MLP tasks | Metadata states all fits complete. Rerun only after fresh baseline parity; preserve pre-frozen markers, seeds and four patient pairs. |
| 6 | Public reviewer DG datasets | 69 original units + colorectal CCA sensitivity1 | 5,928 terminal configurations including excluded Darmanis96; 5,832 eligible for fresh rerun across69 units when Darmanis is excluded. HCL59 tissues account for4,668 conditions. |
| 6 | Public/PTC comparator continuation | Frozen720 unit-method audit entries; planned marker tasks720 | This is a separate incomplete research program; no current generic wrappers for every supervised/pretrained method. Keep its honest pending/excluded statuses. |

The seven-space UMAP is **direct native-R scaled-HVG → uwot**, whereas the original
anchor is **PCA30 → uwot**. Preserve both names; silently substituting the anchor
would change the comparison. PCA/FA/ICA/Isomap/tSNE and KMeans/GMM in that extension
use explicit Python algorithms downstream of identical R matrices. They are
intentional ablations, not the old simplified Python DG-scRNA pipeline. The ICA
truth-blind adaptive convergence repair is frozen in the v8 execution source.

## Exact source/input locations

- Core counts: `data_bench/GSE274546/mtx/<sample>/{matrix.mtx,genes.tsv,barcodes.tsv}`.
  No `obs.csv`, truth or Lfine file is passed to the fit. Hashes for all726 planned
  tasks are recorded in `rerun_inventory.core_tasks.json` from the old input audit.
- Core native inputs/checkpoints: `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/{inputs,GBM}/`.
- Core protocol/order/folds: same root `protocol/{core_tasks.json,sample_order.csv,cohort.csv,patient_folds.csv}`.
- Core marker JSON: same root `markers/libraries.json`; native IDs and all16sets
  preserved. Selected arm `L00_mean` = CM2_glioma_other/mean.
- Native core scientific bodies: `handoff/paper_claim_validation_20260917/{prepare_R.R,score_R.R,terminal.py,legacy_refine.py}`.
  Fixed example derivation: `DGscRNA/research/reference_examples/SOURCE_DERIVATION.json`.
- GBM controls: same code root `geometry_only.py`, `dl_controls.py`,
  `representation_controls.py` plus reviewer `controls/`, `no_clustering/`,
  `embedding/`. Active A1 immutable bundle is under reviewer result
  `source_snapshots/embedding_v8/`; read its manifest before adopting code.
- Lfine provider: `handoff/lfine_compact_20260920/evaluate_lfine.py`, pinned v5
  helper/mapping hashes in `results/hvg_ptc_20260916_v1/lfine_compact_20260920/manifest.json`.
  Current compact proof covers4,719unique conditions, **not all139,392 core arms**.
  A full rerun evaluation must extend the same rule to every arm; do not call this
  smaller view complete grid evaluation.
- Public units/markers/input semantics: reviewer result
  `comparison/analysis_units.csv` and old R campaign `inputs/`, `markers/`.
  `rerun_inventory.json` expands exact69unit paths and marker files.
- PTC grid: `handoff/r_reference_campaign_20260917/ptc_ablation_conditions.json`,
  `ptc_ablation_prepare.R`, `ptc_geometry_only.R`, `reference_score.R`, `terminal.py`.
  Original marker JSON: R campaign result `markers/PTC_original17.json`.
- PTC followups: `handoff/paper_claim_validation_20260917/PTC_FOLLOWUP_PROTOCOL.md`
  and result `PTC_followups/{selection,configurations}/`; final frozen contexts
  must be carried over, not selected again from held-out outcomes.

## Evaluation and provenance contract

- Fit accepts counts, biologically selected marker libraries and sample/batch
  metadata only. Labels are opened by a separate evaluator after outputs freeze.
- Latest GBM report uses compatible-target-set Lfine concordance: frozen native
  panel semantics; present classes with support>=20 excluding Other/nan; all
  cells remain in TP/FP/FN. Do not infer a fine call using a cell's true label.
  Invalid terminal outputs remain missing/invalid, not marker-only results.
- Keep all121/59 and primary97/55 cohorts separate. Candidate choice uses only
  training patients in frozen folds; evaluate on held-out patients. Report valid
  no-op/untrainable terminal states separately from model training.
- Recompute the displayed all-gene UMAP from the actual all-gene native condition;
  old plot_unit.py used an HVG2000 display canvas even for all-gene fits.
- PTC S2 and S3 remain separate original endpoints; TCR detection is not complete
  negative truth. The four patient pairs must stay intact across both groups.
- Xin is RPKM; immune_ALL mixes quantification protocols and uses study-heldout
  groups10X/Freytag/Oetjen/Sun/Villani, not ten independent patients. HCL donors
  keep global fold IDs across tissues. HCL/Baron overlap published scDeepSort
  training sources and are not independent pretrained tests.

## Resource and gate handling

All matrices, fitting, evaluation, plots and numerical comparisons run on SLURM.
Initial profiles are copied from previously successful task families, then checked
against actual new MaxRSS rather than treated as a guarantee:

| Family | Existing profile / caution |
|---|---|
| GBM core | 4 CPU,40GB,4h; infrastructure retry80GB,8h. Start small while existing192-way A1 wave is active. |
| GBM small parity | 4 CPU,32GB; model threads and BLAS bounded. |
| GBM A1 | 4 requested CPU,32GB,24h; nextgen memory coupling may allocate9CPU. Existing6×32waves already use192parallel tasks. |
| PTC CCA preparation | 8 CPU,192GB,18h; original R public preparation128GB,12h. Full-gene jobs may need measured higher profiles. |
| PTC scoring | Fixed-gene64GB minimum: prior32GB jobs OOMed. Full-gene DEG observed~70GB;96/128GB safer starting profile, worker count fixed and monitored. |
| PTC DL | Prior all-gene fits peaked<10GiB; actual jobs retained32GB after rejected scheduler downsizing. Begin32GB,4CPU; do not report proposed16GB as actual usage. |
| New public/PTC method pilots | Public32/64GB4CPU; PTC input192GB8CPU and marker128GB4CPU; these continuation estimates were not resource-tested. |

Nextgen120CPU/node with4027MiB/CPU caps a valid large allocation around448GiB;
do not blindly request480GiB. Check all-user queued-task bounds and current state
before each array release. `rerun_inventory.json` contains the point-in-time queue
snapshot only; refresh it before launch.

Older dispatchers couple acceptance to OneDrive receipts. **Do not invoke their
delivery/finalizer paths.** Latest policy is local notebook only. The new campaign
needs an equivalent hash-verified local GBM completion gate, followed by a PTC
baseline gate; disabling remote delivery does not waive scientific checks.

No package launcher, comparator wrappers or acceptance gate is claimed to exist
by this inventory. Those are implementation tasks for the parent.
