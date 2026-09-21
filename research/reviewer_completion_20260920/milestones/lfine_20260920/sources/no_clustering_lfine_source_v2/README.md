# A2 Lfine re-evaluation v1 — prepared for root review

Read existing GBM seed/terminal predictions only. No normalization, matrices,
clustering, seed scoring or DL refitting. The original A2 fits, labels, L1
metrics/selection/figures and protocols remain archived unchanged.

The endpoint revision follows the latest user request: Lfine only in the compact
notebook. It pins the exact compact provider, legacy v5 helper functions, semantic
mapping, reference labels and existing compact endpoint proof. The broad native
panel → compatible Lfine target-set rule is not a strict one-to-one fine-type
classification. Every cell remains in confusion counts; class averages use the
original support≥20/exclude Other/nan rule. Unknown/Undecided/unmappable calls
receive no compatible-class credit.

Scope: 121 samples × HVG2000/5000 × five saved λ conditions. These are later HVG
sensitivity results, not all-gene main-figure conditions. Fixedλ1 and λ selected
on training-patient mean terminal090 Lfine score are separate. Original patient
folds, exact-tie nearest1 then lowerλ, patient bootstrap2000/seed42, WilcoxonPratt,
and Holm16 contrasts/cohort remain unchanged. All32 contrasts against four
same-budget original routes are emitted. No final winning condition is selected
using held-out outcomes. Initial calls remain diagnostics; terminal070 remains a
sensitivity endpoint. No-op states remain valid; invalid/missing/no-class values
remain explicit and prevent a complete-cohort ranking.

Root first reviews the code and frozen metadata, then writes APPROVED.json with
approved=true plus the exact protocol_sha256 and source_manifest_sha256,
allowed_samples (TKU4163 for the pilot; full roster after review), and
allow_aggregate=false until the full evaluation phase. The
numerical entrypoints refuse to run without this matching local campaign gate.

Pilot command (SLURM only,1CPU4GB sufficient for saved-prediction evaluation):

    python SOURCE/evaluate.py --sample TKU4163

Each sample yields30 A2 rows +24 original-route rows, explicit per-class
TP/FP/FN, and eight original terminal090 comparisons against the existing compact
scores. Every compact field and prediction hash must match, tolerance1e-12.

After root reviews the pilot, full evaluation can use the same command as a
121-task SLURM array (sample names from protocol.samples), or `--all` under an
allocation. `python SOURCE/aggregate.py` requires every sample and blocks silent
sample/patient exclusions. It creates a new summary with all3630 A2 candidate
stage rows,2904 original-route stage rows,20 training-fold choices, fixed/selected
sample results and all32 paired contrasts. Source/inputs/outputs are hash-bound.
No plotting or compact notebook writes occur in these scripts.
