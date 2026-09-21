# Adaptive ICA presentation v1

This is a separate presentation layer over frozen embedding_v8 and unchanged
historical v5/v6 successes. It does not fit models or recompute metrics. Every
ICA entry is named `ICA2 adaptive`; each unit records its actual solver, iteration
cap, and original parallel5000 success/failure. Original figures and summaries
remain in place.

After the normal run/evaluation succeeds, run:

    python render_adaptive_figures.py --tasks /absolute/tasks.json

The helper selects `SLURM_ARRAY_TASK_ID` exactly like run.py; single `space` and
`spaces` tasks are both supported. Non-ICA tasks emit an explicit not-applicable
message and exit successfully. Alternatively pass one existing unit directory:

    python render_adaptive_figures.py /absolute/embedding/SAMPLE/hvg2000/ICA2

All invocations require SLURM. The output is
`UNIT/figures_adaptive_v1/{13 candidate PNGs,13 PDFs,manifest.json,COMPLETE}`.
`COMPLETE` is SHA256 of manifest.json. Manifest `files` maps all 26 names to
SHA256. Other fields bind the fit/evaluation/representation manifests, numerical
policy, renderer source and presentation source manifest. It stores
`space_display="ICA2 adaptive"`, actual_solver, actual_iteration_cap,
canonical_parallel5000_status, fallback_used, and all 13 candidate conditions.
An atomic sibling render-lock directory prevents simultaneous writers. A stale
lock requires review; it is not silently removed.

`verify_manifest(unit_directory)` checks completion, all 26 figure hashes,
source/input/policy hashes, labels and candidate coverage. A dispatcher may use
this function, or implement equivalent checks, for ICA completion. Existing
`figures/` remains unchanged and continues to satisfy the historical output gate.

`select_and_summarize.py` retains the previous primary-cohort boolean assertions,
all consumed evaluation/anchor hash checks and frozen patient-fold selection and
statistics. It additionally requires all ICA adaptive presentations, verifies
old figure hashes, and emits `embedding/summary_adaptive_v1/`. Its candidate,
selection and contrast tables add `space_display`. `ICA_solver_execution.csv`
contains the actual solver/cap/attempt provenance for all 242 ICA units. The
original `summary/` is untouched. All 1,694 units must be complete before any
final ranking is published; no partial ranking or cohort exclusions.
