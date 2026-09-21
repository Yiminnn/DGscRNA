# Final existing-section refresh

Run `run.py --campaign EXTENSION_GATE_2617.json --aggregate CAMPAIGN/evaluation/extensions --selected CAMPAIGN/evaluation/embedding_selected --update-notebook` with the pinned A1 Python inside SLURM. It waits for neither jobs nor partial results: absent or mismatched 726-core, 2,617-extension, aggregate, and independently selected A1 receipts cause rejection before numerical imports. The watcher supplies this command after those stages finish.

The helper updates four existing notebook cells: robustness, neighborhood/training, cellwise seeding, and method comparison. It preserves the other 16 cells, including the completed A1 section and PTC, and preserves core/A1 metadata. The one existing B figure is regenerated; no new figure or chapter is added. It writes `evaluation/final_report/{manifest.json,COMPLETE,notebook_receipt.json,dgscrna_results_candidate.ipynb}` and a pre-update backup. No HTML or remote synchronization occurs.

Scientific scope:

- B: all 408 fresh terminal-threshold rows; 132 unique neighborhood conditions at 0.90, 144 plotted points, 48 descriptive sweeps and 48 seed-aggregated checkpoint rows. Source plotting/reduction code is retained.
- A2: all 2,420 fresh cellwise rows and 1,936 matched original-route rows. Original per-unit CSV parsing, input order, patient means, lambda tie policy and patient folds are retained. All 20 lambda choices and every float64 patient inference input must match exactly before the 32 original confidence intervals and p values can be retained. A one-ULP mismatch blocks completion; sample-level 1e-12 diagnostic tolerance never authorizes retaining different inferential inputs.
- C: 4,719 fresh DG candidate rows from the original HVG2000 partition, 39 candidates per sample. Original round-trip parsing, sorted sample/patient reductions and five folds are retained. All 10 DG fold choices and 114 held-out patient vectors must match exactly. Five non-DG comparator predictions and their conditional inference remain archived, explicitly identified as such. The broader equal24 analysis remains a separate archived record.

A2/C use the accepted per-unit metric CSVs, also bound by completed campaign receipts, to preserve each original parser contract. The generic new aggregate is required and hash-bound. Inference is retained only after exact vector validation; it is not re-estimated here, and these checks do not establish uniform pipeline optimality.

Validation:

- Archived-input function regression: complete original A2/C reduction paths, both one-ULP tamper rejections, B condition mapping/plotting, and four-cell notebook fixture preservation.
- Real newly fitted TKU4163/HVG2000 pilots: 10 A2 candidates, eight original-route anchors and 39 C candidates have exact original metric floats under the specified parsers.
- The real unfinished full campaign is rejected by this helper and the generic aggregator, without data-module imports, output creation or notebook changes.

These tests do not claim that the new full campaign or its report is complete. Full numerical revalidation runs only after its actual completion receipts exist; a mismatch preserves the failed attempt for review.
