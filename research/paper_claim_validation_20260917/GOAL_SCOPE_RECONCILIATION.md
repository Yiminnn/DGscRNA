# Historical recovery task and the current execution goal

The user subsequently clarified: “现在的结果不用完全复线，能差不多就行”.
NEXT_STEPS.md therefore accepts approximate original-workflow reproduction while
requiring the historical discrepancies and metric provenance limits to remain
visible. This explicit clarification supersedes the earlier exact-retraining
stop condition; it does not establish exact historical reproduction.

The original recovery work was actually performed and delivered: native-R
checkpoint scoring and both terminal DL branches ran through SLURM; the complete
archived DL function agrees with the reconstruction at fixed current inputs and
seed42; home/work final labels match Sup for all92,404 cells; all eight original
DG Table2 F1/AUC values were recovered. Per-cell differences were exported rather
than changed to manufacture agreement. Rerun overall F1=.951692/AUC=.938690 is
close to the historical .9518588/.939177. The16NMT and1641TTU label differences,
unavailable historical weights and unresolved V16 Accuracy row remain limitations.

The correction and diagnosis were added to the original notebook and verified in
the authorized original OneDrive directory. Existing receipts:

- results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916/PTC_baseline_correction_upload_receipt_final.json
- results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916/PTC_recovery_diagnosis_upload_receipt.json
- results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/summary/DELIVERY_RECEIPT.json

Closure of the historical recovery task uses the user's later acceptance scope.
It does not mark any pending GBM/PTC experiment, upload, or full-campaign
deliverable complete. The current goal remains execution of NEXT_STEPS.md through
verified GBM delivery, then gated PTC controls, the original notebook, English
workflow figure and verified update of the same OneDrive folder. Preserve all
counterexamples and avoid claiming a predetermined optimal workflow.
