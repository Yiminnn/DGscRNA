# Packaged MLP controls

This research adapter repeats the archived 54 MLP groups: three prespecified samples × HVG2000/HVG5000 × nine configurations. Each group keeps four R routes and two marker libraries at mean cutoff fixed, then produces eight terminal DL conditions and sixteen Lfine threshold rows.

`mlp_controls_tasks.json` freezes the roster. Configurations are the original model, model seeds 0–3 with split seed 42 fixed, widths 128/64 and 512/256, and 5/20 epochs. These are MLP controls; they are not whole-pipeline seed replications.

`refine_controls_packaged.py` is the archived parameterized refinement with only its common-module import replaced. The adapter verifies this exact source transformation and unchanged scientific ASTs, imports the installed package terminal backend, and reads only completed fresh core/pilot scoring outputs. Every extension owns a fresh cache inside its own output directory. Existing extension outputs are refused. Core inputs are checked before and after execution.

All terminal NPZ fields, probability matrices, train/validation splits, model tensors, training histories and valid terminal states must exactly match the corresponding archived `GBM_DL_controls` condition. Marker arms are matched by library and cutoff. Archive models never enter fitting. Lfine evaluation runs only after the exact parity receipt; final completion requires 8 conditions and 16 valid threshold rows. No L1 metrics are produced.

The bounded validation uses the frozen TKU4163/HVG2000 full pilot: task0(original) and task1(model_seed0), each into a separate fresh `extensions/MLP_pilots/` directory. The root authorizes full scheduling only after both tests and the GBM core gate. `mlp_controls_array.sbatch` requests 4 CPUs, 16 GB and 1 hour per group; root controls queue limits and dependencies. No jobs are submitted by the Python adapter.

The public package wheel and the accepted core launcher gate remain unchanged. This adapter and its independent source lock belong to the research extensions, not the published workflow API.

Both bounded pilots passed on SLURM 7204040.274/.273. Each matched all 8 archived terminal conditions exactly (6 trained, 2 valid no-known-label endpoints) and produced 16 valid Lfine threshold rows. Acceptance: `results/hvg_ptc_20260916_v1/package_release_20260920/mlp_controls_pilot_acceptance.json`. The 54-group array remains unsubmitted by this adapter.
