# Reviewer B controls through the released backend

This research adapter preserves the frozen `controls/PROTOCOL.md` scientific
conditions. It does not alter the released package or accepted core outputs.
GBM core remains the prerequisite for campaign units; the explicit TKU4163 /
HVG2000 pilot can use the full-roster package pilot named in `RELEASE_GATE.json`.

## Scope

- **Neighbors:** 66 unique conditions, each with the two original mean-cutoff
  marker libraries. PCA30/UMAP2 SNN k = 10, 20, 40; UMAP neighbors = 15, 30, 60
  with SNN and R HDBSCAN. Duplicate defaults are counted once. All R preparation,
  DEG, density seeds and terminal DL are recomputed, including defaults.
- **Learning:** 18 primary trajectories: TKU4163, NL022, SN040 × HVG2000/HVG5000
  × model seeds 0, 1, 42. Split seed stays 42. The primary marker remains
  CM2_glioma_other/mean; CM2_primary_all_context/mean runs only after a legal
  primary no-training state. Train for 30 epochs and preserve checkpoints at
  5, 10, 20, 30. No truth-based checkpoint or marker selection occurs.

`SOURCE_DERIVATION.json` records the original R and learning sources. Neighbor
R numerical bodies are byte-identical; changes provide explicit input/output
paths and isolated R libraries. Learning model, initialization and training
function ASTs are identical. Its utility import points to the installed package.
The 30-epoch instrumentation remains a research adapter, not a new public API.

## Evidence required

Neighbor acceptance requires exact cell/cluster tables, all 48 initial arms,
all DEG and named density objects, then both terminal arrays, probabilities,
splits, training histories and model weights against archived independent
results. Each unit uses a fresh extension cache and checks that the accepted
source directory is unchanged.

Learning acceptance requires exact per-epoch training/validation history,
checkpoint and final arrays, and model weights against the original learning
experiment. Epoch 10 must also equal the independent original-width same-seed
DL control. Validation targets are held-out **marker pseudo-labels**, not
independent biological truth. Legal no-training states have no invented curve.

The frozen shared evaluator produces Lfine results only after numerical parity:
four threshold rows per neighbor unit and two per saved learning checkpoint.
Unknown remains in the evaluation denominator. Outputs include hash-bound
parity, evaluation and acceptance receipts. Failures preserve their artifacts.

## Execution

Use the gate's installed Python with `-s`, unset `PYTHONPATH`, and allocate four
SLURM CPUs. All matrix reads, model fitting, comparisons and evaluation require
SLURM. The explicit pilots are:

```text
neighbors.py --gate RELEASE_GATE.json --pilot-full-run <gate full TKU pilot>
learning.py  --gate RELEASE_GATE.json --pilot-full-run <gate full TKU pilot>
```

The learning pilot uses seed 42 and 10 epochs; its checkpoints are 5 and 10.
Campaign workers use `neighbors.py --index 0..65` and `learning.py --index 0..17`
with the same `--gate`. They require their corresponding completed pilot and
the newly accepted packaged core for the requested sample/budget. The parent
controller owns all campaign submissions; these scripts do not submit jobs.

Only the explicit TKU defaults are pilot evidence. They do not certify all 66
neighbor conditions or all 18 extended trajectories as rerun.
