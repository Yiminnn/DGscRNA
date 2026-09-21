# Installed-package GBM rerun sources

The published reference package is `dgscrna==2.0.0rc1`:
https://github.com/Yiminnn/DGscRNA/releases/tag/v2.0.0rc1

This archive records the SLURM orchestration, independent verification, frozen
Lfine evaluation and research ablation adapters used to rerun the existing GBM
experiments through that installed release. It preserves the exact execution
sources, including recorded local paths. For a new user's analysis, use the
portable package API and `docs/reference_workflow.md` in the repository.

Scope: 726 core sample/feature units (139,392 terminal configurations), plus 2,617
geometry, MLP, representation, neighborhood, learning, cellwise-seed and seven-space
comparison units. All new terminal results use the original Python DL/refinement
after the original R scientific stages. Marker-only calls and partition metrics
are not substituted for final annotation results.

The source archive and passing pilots do **not** establish full-cohort completion
or optimality. Execution receipts under the local result directory distinguish
pending, failed, valid no-training and trained conditions. Full GBM completion
requires separate core, extension, evaluation and notebook receipts. PTC and
public reviewer reruns follow their own baseline/input contracts and are not
certified by this GBM archive.

The active extension scheduler is `manage_v2.py`. Its separately reviewed
`CONTROLLER_V2_LOCK.json` preserves the original scientific gate and workers.
It handles SLURM's nonzero reply when a throttle update succeeds but completed
array members also produce diagnostics, verifies the actual cap, and records
the preserved state migration. The original scheduler source is retained.

The two UMAP experiments differ intentionally: the original anchor uses PCA30
before UMAP; the seven-space comparison uses scaled HVGs directly. Dataset-specific
preprocessing and the original 2,000-HVG default are preserved. All-gene and other
feature budgets are comparison conditions, not a uniform-optimality claim.

`SOURCE_MANIFEST.json` records byte hashes. Expression matrices, trained weights,
full annotations and notebook binaries are not included in this source snapshot.
