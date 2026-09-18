# Native-R PTC/GBM publication-claim validation

Execution is in progress. This source archive is not a claim that the full experiment or delivery has finished. Read checked output manifests and SLURM accounting for completion.

The accepted scope is in `NEXT_STEPS.md`, with original context in `PLAN_ZH.md` and the exclusion of new Darmanis fitting in `DARMANIS_SCOPE_UPDATE_ZH.md`. `COMPARATOR_PROTOCOL.md` separates matched-marker, broader tuning and labelled-reference comparisons.

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

The scripts contain original site paths under `/fs/scratch/PCON0080/yimin/dgscrna`. This is an execution archive, not a portable Python package release. Raw data, per-cell labels, model weights, installed environments, credentials and rendered notebooks are not copied into Git.

Input audits, source snapshots, package provenance, failed attempts, parity checks and exact-input DL cache signatures are saved under `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/`. The canonical notebook remains `notebooks/dgscrna_results.ipynb`; old cells/results are preserved. No webpage or replacement notebook is created.

## Third-party methods

scType derives its linear score aggregation from the official GPL-3 function; the original source and license are retained. scCATCH3.2.2 and SCINA1.2.0 use isolated installed copies with pinned source archives. SCINA's explicit numerical boundary guard is separately audited and records when the original solver failed; it is not silently described as an unmodified result. SingleR2.8.0 uses labelled training patients under a separate information condition. Pretrained scDeepSort uses its published human-Brain checkpoint when its isolated environment passes validation.

Do not claim a universal optimum. The final report must distinguish an original-default win, an optimized-configuration win, uncertain near-ties and genuine group-specific exceptions.
