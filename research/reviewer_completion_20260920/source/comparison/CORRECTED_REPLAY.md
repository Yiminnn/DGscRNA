# Corrected scDeepSort GBM replay

Scope: frozen 121-sample GBM cohort, published human Brain checkpoint, native
`unsure_rate=2`. This is prediction with pretrained weights, not retraining.
The historical raw-count comparator and its derived evidence remain unchanged.

## Input contract

The official documentation requires Seurat LogNormalize before prediction.
`scdeepsort_corrected.py` starts from the original eligible RNA counts, calculates
`log1p(count / full eligible RNA cell total * 10000)` once, and only then intersects
with the checkpoint gene vocabulary. It preserves floating-point input. Each
sample is checked over all cells and 2,000 native R features. Three pilots also
check every value actually supplied to the predictor against the original R RNA
assay. Normal Brain vocabulary lacks malignant classes; those truth cells remain
in the all-cell denominator. Generic neuronal calls remain ambiguous under the
frozen mapping rather than being resolved using truth.

## Frozen jobs and execution

- Pilot prediction: 7445743 (TKU4163, NL022), 7445852 (largest SN040).
- Independent full-input parity/plot checks: 7445842 and 7445943.
- Remaining 118 samples: array 7446013, initially 2 concurrent, raised to 8 with
  root authorization after three pilots and five cohort samples passed, then
  to 16 after sixteen samples completed without failure/OOM.
- Versioned cohort rebuilding: 7446014, dependent on the full array succeeding.
- Every prediction task requests 4 CPUs, 32 GB, 30 minutes; thread libraries use
  one thread. OSC assigns 9 CPUs to this memory request. No PTC/public jobs.

`scdeepsort_single.py` updates the checksum-verified completion count under a file
lock. `status.json` separates parent array/summary IDs from per-task job IDs.
Task failure is retained and repaired; it never silently removes a sample.

## Versioned outputs

Root: `results/hvg_ptc_20260916_v1/reviewer_completion_20260920/comparison/`.

- `scDeepSort_LogNormalize/<sample>/`: input/prediction manifests, unchanged
  checkpoint/source hashes, native predictions, historical-versus-corrected
  metrics, and exact supplied input CSV. Input CSVs are audit intermediates and
  are much larger than the review figures/tables.
- `corrected_reference_view/`: explicit local symlink view onto frozen unchanged
  inputs and comparator data; new scDeepSort adapters/evaluation and new
  `comparison_summary/`. This is not another raw dataset.
- `evidence_corrected_scDeepSort/`: regenerated six-tool agreement/FMI,
  Unknown/support denominators, patient aggregation and author-L1 context;
  selected-marker panels regenerated under the same frozen expression contract.
- `input_repair_*` and `GBM_scDeepSort_input_repair.*`: descriptive patient-paired
  raw-input versus corrected-input audit.
- `corrected_121_source_verification.csv`, `corrected_comparison_manifest.json`:
  full-cohort source validation, old five-method table parity, independent pilot
  checks and new evidence accounting.

The C inventory (`frozen_dataset_method_matrix.csv` and
`frozen_unit_method_matrix.csv`) is an immutable audit snapshot. Its historical
scDeepSort input finding is resolved by the explicitly separate corrected replay;
other missing comparator conditions remain open. HCL/Baron pretrained weights
have documented training-source overlap and cannot establish independent
external validation on those datasets.
