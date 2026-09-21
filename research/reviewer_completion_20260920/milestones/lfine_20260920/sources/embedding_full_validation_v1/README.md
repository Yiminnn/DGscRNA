# Independent A1 acceptance

Prepared source and metadata only. No scientific validator pilot or full acceptance has run in this version. The frozen original L1 protocol is archival evidence; the canonical notebook uses a separate Lfine endpoint. This package never writes the canonical notebook, producer outputs, dispatcher state, or OneDrive.

The scope is every 121 samples × 2 HVG budgets × 7 representations = **1,694 units**, with **22,022 partitions**, **66,066 initial / terminal 0.90 / terminal 0.70 metric rows**, and all cells retained. It does not represent an all-gene fit. Unit acceptance requires all 13 candidates. Missing or invalid results fail explicitly; they are not dropped or replaced with marker-only output.

## What is checked

- `pilot_baseline.py` retains the independently accepted three-pilot verifier exactly, with only its new isolated output directory changed. This verifies geometry/expression/DL hashes, ordered cells, all partition/score/terminal chains, manual metrics, and original figure hashes.
- `validate_unit.py` adds all per-class and confusion counts, coarse diagnostics, independent contingency-formula ARI/NMI/FMI, terminal NPZ/CSV agreement, known/pool classes, exact recorded split reconstruction, rounded confidence and both thresholds, no-op/trained state rules, parameters, training history, and source whitelist. It never constructs or runs a model. Model files are hashed only.
- Every original PNG/PDF and ICA adaptive figure is checked by manifest, hash and file format. This is an artifact check, not a claim that every figure received visual inspection. Representative visual review remains separate.
- ICA checks record the actual solver, cap, attempt sequence, canonical parallel-5000 result, preserved failure evidence, and deflation acceptance residuals. The frozen source plus recorded residuals is verified; ICA is not fitted again.
- `verify_summary.py` requires all 1,694 proofs and the complete adaptive summary. It independently reconstructs same-budget original-R anchors from final labels, training-patient K choices, held-out sample membership, patient means, 10,000 paired bootstrap intervals, Wilcoxon and Holm adjustment in each 21-contrast family. It verifies the 97-sample/55-patient primary and 121-sample/59-patient sensitivity cohorts, all 420 fold choices and 84 paired contrasts. Inference remains retrospective and conditional on selected predictions.

No scientific module runs on import except `worker.py` and `verify_summary.py`, whose entry points require SLURM. Preparation and sharding use only source files and manifests. All scientific execution commands below must run through `sbatch`.

## Freeze and pilot

The paths are relative to `/fs/scratch/PCON0080/yimin/dgscrna`. The executable is `/fs/scratch/PCON0080/yimin/mamba_envs/dgscrna/bin/python`.

The preparer is one-shot and refuses to overwrite `contract.json`. Once frozen, do not edit package sources. Any repair requires a new version or an explicit pre-execution refreeze; preserve completed proofs.

```bash
python handoff/reviewer_completion_20260920/embedding_full_validation_v1/prepare.py
python handoff/reviewer_completion_20260920/embedding_full_validation_v1/static_check.py
python handoff/reviewer_completion_20260920/embedding_full_validation_v1/make_shards.py --pilot --concurrency 3
```

`make_shards.py` writes a hash-named immutable launch bundle and **prints** an `sbatch` argument list; it never submits. Root reviews and executes that exact command. The three pilot shards validate four already-completed units: hvg2000/PCA2 in TKU4163, NL022 and SN040, plus TKU4163/hvg2000/ICA2. The two TKU4163 units run sequentially in one shard. These pilots exercise the adaptive ICA figure path, verify validator correctness and measure resource use. Do not call an unexecuted script a passed scientific check.

After those jobs pass, root can publish the strict four-proof gate:

```bash
python handoff/reviewer_completion_20260920/embedding_full_validation_v1/accept_pilot.py
```

The gate binds exactly those four current unit proofs and the validation contract. Empty, unrelated or stale proof sets are rejected. Broad sharding and each non-pilot worker require this gate.

## Wave or full validation

```bash
python handoff/reviewer_completion_20260920/embedding_full_validation_v1/make_shards.py --wave-tasks /absolute/path/to/producer_wave_tasks.json --after-job 1234567 --concurrency 8
python handoff/reviewer_completion_20260920/embedding_full_validation_v1/make_shards.py --concurrency 8
```

Use the first form only for an existing frozen producer task list, with its real producing array job ID. The worker verifies the wave-list hash and strict membership in the complete task roster. `afterok` waits for the whole producing wave; it does not manufacture completeness. The second form covers every outstanding unit after production is complete. Running both commands does not submit anything.

Default full release groups at most seven representations per sample/budget, processed sequentially: at most **242 shards**, each **2 CPUs / 8 GB / 4 hours**, suggested combined concurrency **8**. The planner allows at most 16 per bundle; root enforces the combined cap across active bundles and adjusts only after measured pilots. Each unit has an atomic owner lock; duplicate or stale owners fail closed. Retry only after checking the old SLURM job has ended and preserving any failure record; never remove a live lock.

Successful proofs can be reused only when input-manifest, contract and output hashes still match. Scientific worker reuse and the final summary additionally rehash all original data/figure artifacts listed in each proof, so unchanged manifests do not conceal changed terminal predictions. The metadata planner does not perform this scientific acceptance. No existing scientific result is recomputed. Per-unit failures preserve an explicit JSON record; a failure or missing task blocks final acceptance. The checks hash large source binaries but do not load expression matrices; memory use must still be confirmed by the largest SN040 pilot.

After all independent units and the producer adaptive summary finish, root submits:

```bash
sbatch --dependency=afterok:PRODUCER_SUMMARY_JOB:VALIDATION_ARRAY_JOB handoff/reviewer_completion_20260920/embedding_full_validation_v1/summary.sbatch
```

Replace all job placeholders with the real complete dependency set. The summary job is **2 CPUs / 8 GB / 2 hours** and refuses partial input even if a dependency is accidentally omitted. It emits `summary/manifest.json` and `summary/COMPLETE`, not a model or a new producer summary.

## Output contract

All paths below are under `results/hvg_ptc_20260916_v1/reviewer_completion_20260920/embedding_full_validation_v1/`:

- `contract.json`, `tasks.json`, `preparation_verification.json`: pinned scope/source identity and static preparation evidence; not scientific acceptance.
- `units/<sample>/<budget>/<space>/manifest.json` and `COMPLETE`: independently verified unit; metrics, clustering, solver and baseline proof are linked.
- `pilot_baseline/<sample>/<budget>/<space>/validation.json`: unchanged pilot-audit proof in this isolated namespace.
- `PILOT_ACCEPTANCE.json`: exact four-unit gate.
- `summary/manifest.json` and `COMPLETE`: full L1 archival acceptance, including input/proof/source hashes and exact denominators.
- `status.json`: web feed with `work_package=A`, `stage=A1_INDEPENDENT`. Root updates actual job IDs and progress as jobs are submitted; the preparer records no jobs.

## Subsequent Lfine addendum

`lfine_followup_contract.json` is a separate **design-only** interface. It pins the compact v5-compatible evaluator, prefix targets, native-panel semantic mapping audit, and historical parity evidence. The mutable notebook presentation policy is observed only and is not a scientific prerequisite.

After this full archival acceptance, create a separate evaluator and independent verifier with no refitting. Evaluate the same saved initial, terminal 0.90 and terminal 0.70 outputs under the exact compatible-target-set rule: target sets use every observed fine class; macro-F1 averages support ≥20 classes excluding Other/nan; all cells remain in TP/FP/FN, including Unknown and off-vocabulary failures. This is not strict one-to-one fine classification.

**Recompute K from Lfine scores on training patients in the existing folds. Do not reuse L1-selected K as an Lfine optimum.** Apply the new K to the matched initial/0.70 diagnostics; reevaluate the original-R anchor with the same Lfine rule and independently recompute patient statistics. Missing/invalid endpoints remain NA with explicit metric-specific sample/patient denominators; they do not become marker-only results. The original L1 selection and archive remain unchanged. This package has not yet implemented or executed that later Lfine evaluator.
