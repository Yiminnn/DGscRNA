# scDeepSort input and checkpoint audit — 2026-09-20

The current GBM driver `handoff/paper_claim_validation_20260917/deepsort_predict.py`
reads frozen eligible RNA counts, intersects model genes and calls `astype(np.int32)`
before CSV export. Its saved manifests describe raw integer-count input. The
official input contract requires Seurat's default LogNormalize before prediction.
https://scdeepsort.readthedocs.io/en/master/input_requirement.html

The installed v1.0 `deepsort/predict.py` consumes CSV values as graph weights.
It later divides each row by its sum and normalizes edges, but never computes
log1p or Seurat LogNormalize. These internal transforms do not repair raw-count
input. The published paper likewise specifies LogNormalize before the pipeline.
https://pmc.ncbi.nlm.nih.gov/articles/PMC8643674/

Therefore all 121 existing GBM predictions remain historical raw-count runs;
they are not deleted or silently overwritten. The previously completed D/F
supplement correctly describes agreement of those saved outputs, but its GNN
row is not yet a verified correct-input comparison. New corrected predictions
must be separately versioned and then the comparison/consensus supplement
refreshed. Correcting the input does not establish an advantage for any method.

The same scDeepSort paper names HCL as the human training atlas, sourced from
Figshare 7235471, and describes adding Baron GSE84133 to the pancreas reference.
The local HCL and Baron data have these same source identifiers in
`handoff/deck_datasets_provenance.md`. Published checkpoints cannot be counted
as independent external tests on HCL or Baron. Exact barcode overlap has not
been recomputed; the verified source-level overlap already prevents an
independence claim. A donor-excluded retrained GNN would be a distinct design.
https://pmc.ncbi.nlm.nih.gov/articles/PMC8643674/

## Authorized GBM-only corrected-input pilot

First TKU4163 and NL022; no new training. Reconstruct the original RNA assay's
log1p(count / eligible-gene cell total * 10000) on the frozen retained cell set,
check every exported HVG against the native-R DL matrix, and only then intersect
the published model's gene universe. Preserve real-valued normalized inputs;
never int-cast them and never normalize an already normalized matrix again.
Record count/source/model/code hashes, exact cell order, package versions and
gene/label vocabulary. Reuse published Brain weights unchanged and retain the
native unsure_rate=2 threshold. Do not apply DG's 0.90 threshold to another tool.

The Brain checkpoint lacks a malignant class and has generic/fetal neuronal
labels. Keep unsupported classes in the all-cell denominator and do not resolve
generic neurons with author truth. Root may authorize measured full-cohort
execution after these pilots. New PTC/public fitting remains behind the GBM gate.
