# A1 Lfine addendum

Source preparation only; no numerical checks or Lfine runs have been executed in this namespace. This is a separate endpoint for the saved A1 native-R/DL predictions. It never changes fits, native panel labels, the archived L1 results, notebook, OneDrive, or dispatcher.

**The endpoint is the exact compact v5 compatible-target-set Lfine rule.** It credits a native panel's frozen semantic parent against its compatible fine-label set, not a one-to-one predicted subtype. Targets include every observed reference class; macro-F1 averages classes with support ≥20 except Other/nan. All cells remain in TP/FP/FN, including rare classes, Unknown/Undecided and off-vocabulary errors. Every observed class's counts are retained. No reference label chooses a cell's predicted subtype.

The source snapshots preserve `grid_sets.py`, `v5_final_annotations.py` and the compact evaluator. `endpoint.py` AST-extracts only the exact relevant helper functions and uses the frozen CM2_glioma_other parent mapping. The independent checker uses contingency counts and prefix definitions directly, with no call to those metric helpers. The original-R anchor is also checked against the completed compact saved-terminal scores.

## Scope and gates

The frozen scope is 121 samples × hvg2000/hvg5000 × 7 representations × 13 partitions, evaluated at initial, terminal 0.90 and terminal 0.70: 1,694 representations, 22,022 partitions and 66,066 metric rows. This is not an all-gene fit. ICA outputs explicitly retain the adaptive label; actual solver metadata and original accepted proof remain linked.

Root authorized the updated gate on 2026-09-21: each Lfine unit may replay after its **matching original independent unit proof** passes. It does not need to wait for all L1 summaries. The final Lfine selection requires the **complete original archival summary**, all 1,694 original unit proofs, every Lfine metric/count proof, and the complete candidate roster. Invalid or unavailable terminals remain NA; final ranking refuses incomplete candidates. Truth-only lack of a support20 class is recorded separately and cannot be used to choose a favorable patient subset.

K is selected afresh using **Lfine terminal090 training-patient means** in the original saved folds; samples are averaged within patient and training patients weighted equally. Exact ties choose smaller K. The selected Lfine K is applied to matched initial and terminal070 comparisons. L1-selected K is never reused as an Lfine optimum.

Both the 97-sample/55-patient primary and 121-sample/59-patient sensitivity cohorts are preserved. Outputs include per-metric sample counts, total and Lfine-eligible patient counts, all-cell coverage, and truth-derived eligibility. Paired inference uses the same-budget original-R anchor, 10,000 paired patient bootstrap draws with seed20260920, Wilcoxon, and Holm families of 21 contrasts. It remains retrospective inference conditional on selected predictions. `verify_selection.py` independently reconstructs all 420 choices, sample membership, patient means and 84 contrasts.

## Preparation and pilot commands

From the repository root, using `/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python`:

```bash
python handoff/reviewer_completion_20260920/embedding_lfine_v1/prepare.py
python handoff/reviewer_completion_20260920/embedding_lfine_v1/static_check.py
python handoff/reviewer_completion_20260920/embedding_lfine_v1/make_launch.py --pilot --concurrency 3
```

The one-shot preparer freezes hashes; do not edit sources after publication. The planner prints an exact `sbatch` argument list and never submits. Root reviews and submits it. The four proposed pilot units are TKU4163/NL022/SN040 hvg2000/PCA2 and TKU4163 hvg2000/ICA2, grouped into three shards. Each job first runs seven synthetic endpoint boundary tests, then replay and independent count verification. All numerical operations, even those synthetic tests, run only through SLURM.

After measured pilots and review, publish the four-proof gate:

```bash
python handoff/reviewer_completion_20260920/embedding_lfine_v1/accept_pilot.py
```

For later already accepted units, `make_launch.py --tasks /absolute/path/to/unit_tasks.json --concurrency 8` prints a bounded launch. The JSON list uses exactly `sample`, `budget`, `space`. Use `--after-job JOBID` only with a real dependency; every unit still checks original acceptance at runtime. No planner submits jobs. Each group runs at most seven units sequentially, 2 CPUs / 8 GB / 2 hours; root controls combined concurrency across bundles and adjusts based on actual pilot measurements. Duplicate/stale locks fail closed.

After all full acceptance gates are satisfied, root can submit `summary.sbatch` with the complete dependency set. It runs the selector and then its independent verifier, 2 CPUs / 8 GB / 1 hour. The selector's `COMPLETE` is explicitly **pending independent selection verification**; final acceptance is only `summary_verification/manifest.json` plus its matching `COMPLETE`.

## Outputs

Under `results/hvg_ptc_20260916_v1/reviewer_completion_20260920/embedding_lfine_v1/`:

- `protocol.json`, `contract.json`, `shards.json`, `preparation_verification.json`: source/scope preparation, not scientific acceptance.
- `evaluation/<sample>/<budget>/<space>/`: 39 A1 rows, three anchor rows, all observed-class counts, eligibility and input/output hashes.
- `verification/<sample>/<budget>/<space>/`: independent count and exact compact-anchor checks.
- `summary/`: every candidate, freshly selected Lfine K, selected sample rows, patient means/counts, and paired statistics.
- `summary_verification/`: independent full selection/inference acceptance.
- `status.json`: A1_LFINE web progress; no stage is marked complete before its proof exists.

The mutable notebook presentation policy is not a scientific hash prerequisite. Current display instructions remain Lfine only, actual all-retained-gene main figures, and HVG comparisons later; the separate notebook session controls that delivery.
