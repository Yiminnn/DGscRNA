# GBM/PTC research source snapshot — 2026-09-17

This branch preserves the current GBM experiments, cluster-figure generation, PTC R-workflow reconstruction, original-output reconciliation and DL initialization diagnostics. The package baseline is commit `4bf17c4cb9518ab427fa92e8bf5bb4d594d30e00` from `align-r-reference`.

All additions are confined to this research directory. The Python package, existing examples, tests, CI configuration and other branches retain their previous contents. This is a source archive of the executed research workflow, not a new package release or a claim that the simplified Python package reproduces the paper.

## Navigation

| Location | Purpose |
|---|---|
| [Python/R comparison](PYTHON_R_DIFFERENCES_ZH.md) | Differences, verified code issues, platform effects and alignment priorities |
| [Source manifest](SOURCE_MANIFEST.json) | Original relative paths, SHA256 checksums and local Python dependencies |
| `source/handoff/hvg_ptc_20260916/` | GBM feature/representation/clustering experiments and terminal annotation evaluation |
| `source/handoff/gbm_cluster_figures_20260916/` | Complete GBM cluster and annotation figure atlases |
| `source/handoff/ptc_recovery_20260916/` | Historical-file inventory, R reconstruction and earlier two-group experiments |
| `source/handoff/ptc_paper_baseline_20260916/` | Corrective original-checkpoint R/DL reruns, Sup/Table2 reconciliation and initialization diagnostics |
| `source/handoff/repo_snapshot_20260917/` | Source packaging and synthetic review probes |
| `source/handoff/gbm/`, `source/handoff/g274_table4/`, other included modules | Dependencies used by these research scripts |
| `reference/ptc_archive/` | Byte-preserved R and Python source from the recovered PTC archive |
| `validation/` | Synthetic audit results and validation records; no patient data |

## Scientific state at this snapshot

- GBM experiment and figure delivery is complete. These Python experiments have their own documented settings; they are not equivalent to the complete original R workflow.
- Original PTC uses a common eight-sample CCA checkpoint, then selects Thyroid/Seurat/none/final-DL for NMT and Pubmed34663816/UMAP-HDBSCAN/mean/final-DL for TTU. The earlier independent two-group CCA scripts remain dated experimental code, not the corrected paper baseline.
- Original final Sup labels and all eight DG F1/AUC values have been recovered. Exact historical terminal retraining and the V16 Accuracy source remain unresolved.
- Fixed-input legacy DL reconstruction has been verified against the complete archived function. This verification applies to the reconstructed inputs and fixed seed, not to the simplified package or to unknown historical model states.
- Later PTC ablations remain paused. “RL” in the work request means terminal DL/refinement.

## Running and portability

The copied scripts retain their original bytes, including site-specific workspace, environment and SLURM paths. Their original workspace layout is recorded in the manifest. To inspect the version, read this directory directly. To run it elsewhere, first configure an explicit working copy of those paths, input locations, environments and SLURM account/partition; preserve this source snapshot for comparison. A clean repository checkout alone does not contain the research inputs and is not sufficient to rerun the cohort analyses.

The code expects separately held GBM count matrices, marker libraries and prepared-input manifests, and PTC archive/checkpoint/Sup files at the paths named in the scripts. Generated intermediate files and some third-party dependencies are also external. No expression matrices, patient metadata, per-cell labels, model weights, rendered notebooks or OneDrive contents were added here.

All scientific computation, numerical tests, matrix loading, fitting, plots and large archive extraction in this work use SLURM. The archived orchestration and delivery scripts have external side effects when run; inspect their inputs and destinations before executing them in another workspace. This repository update does not execute cohort jobs or publication scripts automatically.

For a future reusable implementation, establish an explicit R-reference mode and a separately named Python-native mode, with persisted expression layers, feature lists, partitions, initial labels, model state and numerical parity checks. The accompanying comparison describes the work required; this branch does not silently change either algorithm.
