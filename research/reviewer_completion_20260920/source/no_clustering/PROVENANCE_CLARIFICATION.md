# A2 provenance clarification

Added after the frozen pilot, without changing candidate definitions, fitted outputs,
selection, or the frozen scientific source.

“Original” in A2's regression and comparator descriptions identifies the native-R
reference workflow reproduced in `paper_claim_validation_20260917`. It does not
claim byte-identical historical manuscript model weights.

The refinement architecture and operations are unchanged relative to that frozen
reference reproduction. The shared reproduction helper explicitly sets Python,
NumPy, and PyTorch initialization to 42 and sets the known-cell split generator to
42. The pilot regression proves that A2's copied helper reproduces the reference
runner's outputs. It does not prove that the historical manuscript training used
global initialization seed 42; the local manuscript/source inventory established
the split seed but did not establish that historical global initialization.

Accordingly, the frozen protocol's phrase “original ... initialization seed42”
should be read as “the same explicit initialization seed 42 as the reference
reproduction,” with historical initialization unverified. This clarification is
reporting provenance only; every A2 candidate and its original-route comparator
retains the same frozen reproduction settings.

Evidence:

- `results/hvg_ptc_20260916_v1/reviewer_completion_20260920/no_clustering/source_v2/legacy/legacy_refine.py`, lines 62 and 71: global and split seeds.
- `results/hvg_ptc_20260916_v1/reviewer_completion_20260920/no_clustering/regression/manifest.json`: executed regression, SLURM 7445794.
- `handoff/reviewer_completion_20260920/manuscript/SOURCE_GAPS.md`: original-source limitations.
- `results/hvg_ptc_20260916_v1/reviewer_completion_20260920/no_clustering/source_v2/SOURCE_MANIFEST.json`: immutable scientific source fingerprints.

The A2 intervention is a replacement of cluster-DEG seed construction with
cell-wise marker-expression seeding. It is not a pure deletion of clustering,
and its results cannot establish the performance of every possible method that
does not cluster cells.
