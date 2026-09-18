# Native-R PTC/GBM publication-claim validation

Execution is in progress. This source archive is not a claim that the full experiment or delivery has finished. Read checked output manifests and SLURM accounting for completion.

The accepted scope is in `NEXT_STEPS.md`, with original context in `PLAN_ZH.md` and the exclusion of new Darmanis fitting in `DARMANIS_SCOPE_UPDATE_ZH.md`. `COMPARATOR_PROTOCOL.md` separates matched-marker, broader tuning and labelled-reference comparisons.

`PTC_FOLLOWUP_PROTOCOL.md` specifies the bounded continuation: reuse the completed
native-R grid, evaluate the verified four patient pairs, separate initialization
from PCA/UMAP seed effects, and test a uniform marker-retention scoring universe.
Every new PTC scientific entry point checks the verified full-GBM delivery receipt.
Preparing the adapters does not imply that the PTC parity pilots or new fits passed.

## Scientific contract

- GSE274546: original121samples/59patients,429305cells; fixed primary97samples/55patients. Counts and author labels are stored separately. All prediction rows retain exact cell IDs.
- Native Seurat RNA branch: LogNormalize10000, VST500/1000/2000/3000/5000/all, ScaleData, PCA30, uwot UMAP2, SNN/Louvain0.5 or R HDBSCAN minPts50. Single-sample RNA DEG scoring uses all eligible RNA genes; DL uses normalized selected genes. These stage gene lists are exported explicitly.
- Six feature budgets x four clustering routes x16marker libraries x3density cutoffs. Native panel labels are kept through the historical10-epoch MLP; the final result is terminal DL/refinement, including its true trained/no-op/untrainable status.
- Core confidence0.90 and probability-derived0.70 sensitivity; Unknown and unmapped calls count as errors. L1 present-class macro-F1 is primary, fixed11secondary. Clustering ARI/NMI is reported separately.
- Geometry-only controls fix full-RNA scoring and normalizedHVG2000DL inputs. MLP-only controls fix representation, clusters, seed labels and training split; native defaults must reproduce every class probability before variation.
- Patient folds group every sample of a patient. Marker/configuration selection uses training patients; heldout labels only score. The cohort was previously explored, so this is retrospective validation.
- The historical PTC NMT PCA-SNN branch and TTU UMAP-HDBSCAN branch are preserved. No independent Pu/GSE184362 cohort or new Darmanis fit is introduced.

## Execution

All scientific computation and plotting require SLURM. `job.sbatch` sets the site environments; `dispatch.py` schedules core fits, terminal checks, evaluation and all-route figures. `aux_dispatch.py` schedules controls and audited marker comparators. `finalize_launcher.py` freezes the current finalizer only after all726core units pass; `core_finalize.py` updates the existing notebook and verifies the upload to the already authorized OneDrive directory. Further controls have separate completion gates.

`scina_library_dispatch.py` gates each library on its own three size pilots before
parallel cohort release; evaluation still requires all 16 libraries.
`full_GBM_launcher.py` waits for all controls, competitors and measured resource
points before the complete GBM report/delivery. `ptc_dispatch.py` then releases
saved-grid evaluation and native-default parity pilots, followed by 22 controlled
PTC preparation units and 50 MLP task groups. It requires a marker-retention resource
pilot before expanding that stage. Scientific failures require log review; they
are never turned into zero performance or silently dropped.

The scripts contain original site paths under `/fs/scratch/PCON0080/yimin/dgscrna`. This is an execution archive, not a portable Python package release. Raw data, per-cell labels, model weights, installed environments, credentials and rendered notebooks are not copied into Git.

Input audits, source snapshots, package provenance, failed attempts, parity checks and exact-input DL cache signatures are saved under `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/`. The canonical notebook remains `notebooks/dgscrna_results.ipynb`; old cells/results are preserved. No webpage or replacement notebook is created.

## Third-party methods

scType derives its linear score aggregation from the official GPL-3 function; the original source and license are retained. scCATCH3.2.2 and SCINA1.2.0 use isolated installed copies with pinned source archives. SCINA's explicit numerical boundary guard is separately audited and records when the original solver failed; it is not silently described as an unmodified result. SingleR2.8.0 uses labelled training patients under a separate information condition. Pretrained scDeepSort uses its published human-Brain checkpoint when its isolated environment passes validation.

Do not claim a universal optimum. The final report must distinguish an original-default win, an optimized-configuration win, uncertain near-ties and genuine group-specific exceptions.

Reviewer C1 explicitly requests repeated runtime measurements. Resource execution
therefore includes three independent process runs per method and size (45 total),
with identical nested counts, isolated result/DL caches, mean and sample SD, and
actual scheduler allocations. OS file-cache state remains uncontrolled.

`analysis_dispatch.py` releases independent summary/diagnostic jobs as their own
prerequisites finish. `unknown_expression.py` adds patient-paired, within-author-class
Unknown expression contrasts; they are exploratory associations because annotation
uses the same expression. `ptc_comparator_replay.py` reuses corrected SignacX outputs,
separates explicit-T CellStates from coarse TNK, and corrects old endpoint framing
without fitting another model. It is subject to the same GBM-before-PTC gate.
`reviewer_evidence.py` links the PI/reviewer requests to verified artifacts while
retaining manuscript, deposition and biological-validation limitations.

Metadata controllers use `metadata_chain.py`: one-hour allocations renew between
complete ticks after50minutes, leaving fitting jobs untouched. This avoids multi-day
metadata allocations being backfilled days later. The nextgen memory/CPU coupling
is also enforced: MaxCPUsPerNode120 and MaxMemPerCPU4027MiB make480GiB infeasible;
resource submissions cap at448GiB, with measured-pilot headroom and actual allocation
recorded. Scientific failures still require review; a bounded pilot recovery only
retries a scheduler-confirmed timeout/OOM/node failure and preserves the attempt.
