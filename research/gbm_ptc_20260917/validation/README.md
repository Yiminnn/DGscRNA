# Validation

SLURM job `7340233` completed with exit `0:0`. It ran four synthetic probes against the unchanged package at `4bf17c4` and then the existing `tests/test_r_alignment.py` suite (six tests passed).

The probes verified missing preprocessing expression-layer preservation, overwritten DEG state across multiple partitions, unused Harmony coordinates in the HDBSCAN PCA route, and `Unknown` entering the training vocabulary. Two probes intentionally mock algorithms to isolate data routing. These findings are not patient-level effect sizes and do not constitute a failed cohort experiment.

`package_probe_results.json` records the package versions, probe definitions and observed behavior. `snapshot_validation.json` records source hash, Python/R/shell syntax and archive-scope verification. Scientific computation used SLURM; packaging copies and Git operations only handle source files and metadata.

Passing the existing six tests does not imply that the complete Python pipeline matches the original R workflow. This branch records these uncovered gaps without changing the established package or experiment implementations.
