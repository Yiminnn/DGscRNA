# Fixed-HVG2000 DL geometry control

`geometry_control.py` is a research adapter around the installed package terminal backend. It preserves the original 363 units: 121 samples × geometry budgets HVG2000/HVG5000/all; each unit contains four routes × two prespecified libraries (`CM2_glioma_other`, `CM2_primary_all_context`) at the mean seed cutoff.

Geometry-specific initial calls come from the newly accepted core run. DL always receives that sample's normalized RNA HVG2000 expression. The all-RNA scoring gene universe remains fixed; cluster-specific DEGs and initial calls belong to the corresponding geometry condition. The adapter performs no new DEG calculation and no truth-dependent marker selection.

Requirements:

- The actual published release gate and installed runtime must pass `core_task.validate_gate`.
- Normal operation requires accepted new core outputs for the geometry budget and HVG2000, including their exact-parity and Lfine receipts.
- An explicitly named pilot may instead use only the frozen full-roster TKU4163/HVG2000 pilot. This does not count as a completed 363-unit campaign.
- Run through SLURM with four CPUs, using the frozen runtime's Python `-s`, without `PYTHONPATH`.

Pilot command, only after the root has published the gate:

```bash
/absolute/runtime_v3/bin/python -s handoff/package_release_20260920/geometry_control.py \
  --gate /absolute/RELEASE_GATE.json \
  --pilot-full-run /absolute/pilots/v3_TKU4163_hvg2000_full
```

For one authorized campaign unit, use `--index N`, where N is 0–362 in the frozen core task roster filtered to HVG2000/HVG5000/all. The adapter submits no jobs itself.

Outputs go to `<output_root>/extensions/geometry_fixed_DL2000/<sample>/<budget>/`; the explicit pilot goes to `<output_root>/pilots/geometry_fixed_DL2000_TKU4163_hvg2000/`. Existing attempts cause an error and are preserved. Score files are copied into the extension. Only the fixed expression binary may be hardlinked and the backend opens it read-only. Core directory inventories and all shared/copied input hashes are checked after fitting.

Each extension owns a fresh DL cache; no model/cache is copied from core or historical results. Byte-identical inputs may share a fit only within the new extension unit. The independently frozen verifier compares every NPZ array, final prediction, probability, training split/history and model weight exactly against the original `terminal_geometry_only_DL2000` outputs. Archived outputs are opened for verification after fresh fitting. No old L1 evaluation is imported or displayed.

`geometry_parity.json` records exact verification before evaluation. The shared frozen `evaluate_extension_lfine.py` then evaluates only Lfine, producing 16 rows: eight terminal conditions at two thresholds. `GEOMETRY_VERIFIED_COMPLETE` is written only after exact parity, complete Lfine evaluation and unchanged core inputs are verified. `geometry_acceptance.json` binds both proofs and records the source gate, installed package, source-run receipts, fixed expression hash, terminal states, fresh-cache provenance and numerical tolerance. A failure writes `failure_preserved.json`; it never changes parameters or modifies package v3.

Implementation has not submitted a geometry job. Root owns pilot/full-campaign launch authorization.
