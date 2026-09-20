# Source gaps and correction ledger

Scope: local document and metadata inventory, 2026-09-20. No new PTC processing, wet-lab inference or public submission occurred. Paragraph IDs refer to `original_v16.paragraphs.txt`; they are XML positions, not page numbers. The original DOCX and existing response remain unchanged. Full file hashes are in `source_inventory.json`.

| Item | Conflicting or available source | What can be stated now | Remaining source / acceptance |
|---|---|---|---|
| Viability | V16 P0558: >80%; P0561: >85% | The manuscript contradicts itself. Original delivery table documents samples/counts, not an independently verified viability protocol. | Original dissociation SOP, recorded per-sample viability or laboratory confirmation tied to an original record. Do not pick one threshold based on plausibility. |
| Strainer | V16 P0558: 70 µm; P0561: 40 µm | Both appear in Methods. A two-stage filtration procedure is possible but is not documented by this conflict alone. | Original laboratory protocol and whether sizes differed by step or sample. Do not invent sequential filtration. |
| GEX chemistry | V16 P0562: 3′; six original GEX reports: Single Cell 5′ PE | Reports support 5′ PE for MT-2, N-1, TU-1, TU-2, T-1, T-2. Original GEX reports missing for MT-1/M1 and N-2/N4. | Original kit/order/run records and missing reports; do not extrapolate six reports to the remaining two without evidence. Distinguish spatial 3′ chemistry from GEX. |
| TCR pipeline version | V16 P0662: Cell Ranger 7.1; all eight original VDJ report metadata: 3.0.1 | Report metadata contradict the manuscript version. | Preserve original reports and establish whether any later rerun generated the analysis input. If not, revise to the report-supported version; do not infer a rerun. |
| GEX pipeline versions | Five available reports: 3.0.1; T-2 report: cellranger-5.0.0; MT-1/N-2 reports missing | Versions are library-specific, not one verified universal version. | Link each final count matrix to its report and hash; obtain missing report sources. |
| Historical DoubletFinder | V16 P0599 and archived writing note line226: 2.0.3; V16 resource row P0486: 2.0.6; P0487 attributes Xun | Attribution should be reconciled with the original McGinnis reference P0430. The recovery environment is 2.0.6 but was inspected in 2026 and is not historical-run evidence. | Original sessionInfo/lockfile/install log or saved run provenance. The helper API suffix difference alone does not identify a unique historical version. |
| Historical Seurat | V16 P0483: v5.3.0; original saved object inventory: Seurat 4.0.1 | Current replay version and historical saved-object version must be separate. | State the version associated with each archived object and each refit; inspect original run provenance rather than assigning V5 retrospectively. |
| Random initialization | V16 P0666 says all stochastic operations, including initialization, seed42; original helper records split seed42 but no global seed proof | Reproducible new seeds can be recorded; the original split seed does not establish seeded network initialization. | Preserve uncertainty; report absence of original model weights and use accepted refit tolerance. |
| Delivered-to-analysis count | V16 P0563 calls110,497 post-QC; P0038/P0387 and original final endpoint use92,404 | 110,497 is verified delivered GEX barcode count.92,404 is the later preserved historical endpoint. | Join barcode membership through original stages; never label every missing barcode a doublet or QC failure by subtraction. |
| Raw reads and accession | Existing GEO package README and submission_gaps.tsv record no supplied FASTQ/BAM/CRAM and no created accession | A local processed-file preparation package exists. This bounded audit did not search every home/archive location and does not assert raw reads cannot exist elsewhere. | Original raw-read source or documented handling with GEO; actual accession and reviewer-access verification through an authorized account. |
| Historical Accuracy | V16 Table2 vs recovered metric reconciliation | Historical F1/AUC are reconstructable under named original definitions; original Accuracy source/definition remains unresolved. | Original computation script or demonstrable definition. Preserve paper-reported historical values separately; use recomputable definitions for new tables. |

Primary metadata source: `results/hvg_ptc_20260916_v1/geo_submission_v1/PTC_GEO_submission_20260916/00_metadata/pipeline_metadata_from_original_reports.tsv`. It records original archive member paths per library. Original GEX reports are preserved under that package's `01_scRNA/original_received/`; original VDJ reports are in `03_TCR/original_received/`. This draft reads existing extracted metadata; it does not newly validate every report-to-matrix linkage.

The recovery-version observation is from `results/hvg_ptc_20260916_v1/ptc_recovery/inventory_r/environment.json` and `sessionInfo.txt`, explicitly a later environment. The historical object evidence is `results/hvg_ptc_20260916_v1/ptc_recovery/inventory_workspace/objects.json` and `object_01.commands.txt`.

## Cell-flow ledger: known entry counts and missing stage attribution

These are copied from the existing verified `00_metadata/original_delivery_table_reconciliation.tsv`; no new matrix counts or attrition inference were performed.

| Analysis sample | Original GEX ID | Delivered GEX barcodes | QC/singlet/final membership status |
|---|---|---:|---|
| MT-1 | M1 | 7,571 | Per-barcode stage join pending |
| MT-2 | M2 | 19,264 | Per-barcode stage join pending |
| N-1 | N3 | 9,653 | Per-barcode stage join pending |
| N-2 | N4 | 18,077 | Per-barcode stage join pending |
| TU-1 | P1 | 9,936 | Per-barcode stage join pending |
| TU-2 | P2 | 17,097 | Per-barcode stage join pending |
| T-1 | T3 | 19,677 | Per-barcode stage join pending |
| T-2 | T4 | 9,222 | Per-barcode stage join pending |

The documented full delivered total is110,497; the preserved downstream annotation total is92,404. The arithmetic difference is not a filtering diagnosis. Samples map to patients through the original delivery table and `sample_crosswalk.tsv`; sample identity must accompany every barcode join. Some TCR barcode suffixes differ from GEX suffixes. Retain the original strings and document any within-sample nucleotide-core reconciliation rather than silently stripping sample identity.

Candidate stage sources already located:

- Delivered barcode files: GEO package `01_scRNA/original_received/` and their archive member paths in the file manifest.
- Historical endpoint and preserved object metadata: `ptc_recovery/inventory_workspace/objects.json`, `object_01.commands.txt`; S2/S3 reconciliation in `ptc_paper_baseline/home_original_sup_manifest.json` and `original_final_sup_references.csv.gz`.
- Recovered fixed analysis inputs: `ptc_recovery/fixed_QC_inputs/cells_QC_sample_patient.csv`, `MTN.cells.csv`, `TUT.cells.csv`. Their filename does not by itself prove the reason for an excluded barcode.
- Original intermediate candidates: `ptc_recovery/archive/tcr/ptc_val/processed_data/PTC_subt_final_metadata.csv`, `ptc_batch/thyroid_sample_counts.csv`, and archived Seurat checkpoints inventoried in `inventory_workspace/`. A T-cell subset checkpoint cannot substitute for a complete pre/post-QC checkpoint.
- Original helper QC logic: archived `DGscRNA-Share/R/source.R:94-151`. The later `science_22_pipeline.ipynb` loads already annotated data before its later QC and must not be retroactively treated as the QC that created the historical endpoint.

Required ledger fields: `sample_id`, `patient_id`, `barcode_original`, `barcode_join_key`, `delivered_present`, `object_creation_present`, `qc_present`, `singlet_present`, `final_S2_present`, `final_S3_present`, `first_observed_exclusion_stage`, `recorded_reason`, `reason_source`, `source_file_hash`, `status`. An unavailable intermediate checkpoint is an explicit missing-source state, not an inferred reason.

## Existing response draft is a separate historical document

`paper/revision/Review_Reply.docx` is dated2026-07-28 in its extracted P0003. Its P0008 describes eleven external datasets and a scaling curve to584,207 cells, while subsequent sections contain simplified-Python branch results and older reviewer IDs. These statements cannot be carried into the current original-R revision as newly validated evidence. The new response draft follows the approved35-row mapping and links each point to existing or pending evidence, leaving the older reply intact.
