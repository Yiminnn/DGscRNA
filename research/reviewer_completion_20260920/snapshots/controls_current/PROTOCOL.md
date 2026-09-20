# Reviewer B — GBM controls (frozen 2026-09-20)

Execution scope: TKU4163, NL022, SN040 (existing size-selected pilots), with HVG2000 and HVG5000. This is not a full-cohort sensitivity claim. Existing original R objects, all-RNA density scorer, normalized selected RNA DL input, markers, cutoffs, cell order and truth mapping are unchanged.

## Neighbors

66 unique conditions (11 per sample/budget): PCA30 and PCA30→UMAP2 with SNN k10/20/40; PCA30→UMAP2 n.neighbors15/30/60 with SNN and R HDBSCAN. The identical SNN k20 / UMAP n30 condition appears once. All other knobs remain original: SNNresolution0.5, minPts50, PCA30, UMAPcosine/min.dist0.3, geometryseed42. Two predeclared marker contexts: CM2_glioma_other/mean and CM2_primary_all_context/mean. Terminal DL modelseed42, split42, epoch10. Every unique condition has clustering plus terminal annotation figures; the fixed HVG2000 coordinates are explicitly only display coordinates.

Default-parity gate: independent new prepare/score execution for TKU4163/HVG2000 PCA30-SNN, UMAP2-SNN and UMAP2-HDBSCAN. Require identical cluster labels, all48 initial marker/cutoff arms, and both fixed-context terminal labels, class probabilities and training split. Remaining default conditions reuse immutable original outputs only after this gate. Changed neighbor conditions recompute their own R clusters, DEG and density seeds; terminal cache is restricted to byte-identical inputs, labels and parameters in the new controls directory.

## Real validation curves

18 primary learning conditions: three pilots×two budgets×modelseeds0/1/42. The historical split42 and initialization/training sequence stay unchanged. Actual per-epoch training loss/accuracy and held-out marker-pseudo-label validation loss/accuracy are recorded for30epochs. Save models, all-class probabilities, predictions and0.90/0.70 endpoints at5/10/20/30epochs. For every trainable condition, epoch10 class probabilities, labels and train/validation cells must exactly equal the existing original-width same-seed control. A small seed42/10epoch parity run precedes the array.

Primary marker CM2_glioma_other/mean, route UMAP2_HDBSCAN_R. Only a legal no-training primary endpoint permits an additional separately named CM2_primary_all_context/mean condition. Preserve no-op/invalid states; no fabricated curves. Unknown stays in the evaluation denominator. Truth is opened only after predictions freeze. No early stopping, marker selection or checkpoint selection uses author truth. Validation targets remain pseudo-labels, not independent biological labels.

## Artifacts and compute

New source only under this directory; old sources remain read-only. Results: results/hvg_ptc_20260916_v1/reviewer_completion_20260920/controls/. All computation and figures through SLURM nextgen/ascend-default, accountPCON0080,4CPU/32GB. Array concurrency≤2. Atomic status.json is refreshed after each condition for the progress webpage. Parity and all task checks must pass before summary/COMPLETE. PTC remains outside this bounded package until the parent's GBM gate.

2026-09-20 execution update: after parent verified current project/QOS headroom and pilot resource use, array concurrency was raised from2to6 using ArrayTaskThrottle. Each fit remains4CPU/32GB with identical scientific parameters and the same worker count. This supersedes only the initial≤2 launch throttle.
