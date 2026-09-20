# Original supplements: inventory and remaining verification

Inventory date: 2026-09-20. Source files are preserved in `paper/submission_v16/`.
Document extraction and container checks ran in SLURM7204040.43. XLSX CRC and
workbook-XML checks establish container readability, not correct biological
content or final spreadsheet-reader rendering. No matrix loading or new PTC
computation occurred in this documentation task.

The original V16 supplement descriptions are in extracted paragraphs P0385-P0389.
New Supplement A/B must supplement these requested original files, not silently
replace their identity or endpoints.

## Supplementary Figure S1.pdf

Original path: `paper/submission_v16/Supplementary Figure S1.pdf`

SHA-256: `485541fd83d848270878a57c3438331287de7cf505c942e6076c7890724af7b1`

Bytes: 786655.

Remaining: Final figure legibility/caption/numbering and reader rendering need verification.

## Supplementary Figure S1.png

Original path: `paper/submission_v16/Supplementary Figure S1.png`

SHA-256: `7e4223e4a03df5ea8a4c7e50af88ee7b04977709c7776cbc07a5496ab99cef4d`

Bytes: 4405037.

Remaining: Final figure legibility/caption/numbering and reader rendering need verification.

## Supplementary Figure S1.tif

Original path: `paper/submission_v16/Supplementary Figure S1.tif`

SHA-256: `30d1d8e42a7e5cd21830d5f3a74cef44ecbf58145d2025d0cf9bef6738804ecb`

Bytes: 6310052.

Remaining: Final figure legibility/caption/numbering and reader rendering need verification.

## Supplementary table S1_Cell_markers.xlsx

Original path: `paper/submission_v16/Supplementary table S1_Cell_markers.xlsx`

SHA-256: `00db62499fa71f83013136b243ead775faeaf51340c24ceb76c69d7e41146aeb`

Bytes: 94196.

Container CRC: passed. Workbook sheets: `Supplementary Table 2`; `CellMarker2.0_Lymph`; `CellMarker2.0_Lymph node`; `CellMarker2.0_Lymphoid tissue`; `CellMarker2.0_Thyroid`; `CellMarker2.0_Epithelium`; `CellMarker2.0_Blood`; `CellMarker2.0_Thymus`; `CellMarker2.0_AllTissues`; `NCOMREFF`; `HPA_CT enhanced_Protein`; `HPA_CT enriched_Protein`; `HPA_Group enriched_Protein`; `HPA_CT enriched_Other`; `HPA_CT enhanced_Transcript`; `HPA_CT enriched_Transcript`; `HPA_Group enriched_Transcript`; `HPA_AllCategories`.

Declared XML sheet dimensions: `xl/worksheets/sheet1.xml: A1:A2`; `xl/worksheets/sheet2.xml: A1:K10`; `xl/worksheets/sheet3.xml: A1:G18`; `xl/worksheets/sheet4.xml: A1:I14`; `xl/worksheets/sheet5.xml: A1:M11`; `xl/worksheets/sheet6.xml: A1:R18`; `xl/worksheets/sheet7.xml: A1:BF311`; `xl/worksheets/sheet8.xml: A1:Q61`; `xl/worksheets/sheet9.xml: A1:BF437`; `xl/worksheets/sheet10.xml: A1:I9`; `xl/worksheets/sheet11.xml: A1:CY14`; `xl/worksheets/sheet12.xml: A1:G7`; `xl/worksheets/sheet13.xml: A1:AH12`; `xl/worksheets/sheet14.xml: A1:B3`; `xl/worksheets/sheet15.xml: A1:BD14`; `xl/worksheets/sheet16.xml: A1:ED6`; `xl/worksheets/sheet17.xml: A1:AE11`; `xl/worksheets/sheet18.xml: A1:IZ14`. These are metadata, not newly computed cell counts.

Remaining: Full spreadsheet reader, endpoint/content validation and submission cross-reference check remain; XML dimensions are metadata, not recomputed cell counts.

## Supplementary table S2_T cell type identification_DGscRNA.xlsx

Original path: `paper/submission_v16/Supplementary table S2_T cell type identification_DGscRNA.xlsx`

SHA-256: `273c57a6234545a889ab2221fd43fc9ae9c077c76d1d25b81370132e0e670fe8`

Bytes: 10407244.

Container CRC: passed. Workbook sheets: `Supplementary Table 2`.

Declared XML sheet dimensions: `xl/worksheets/sheet1.xml: A1:Q92408`. These are metadata, not newly computed cell counts.

Remaining: Full spreadsheet reader, endpoint/content validation and submission cross-reference check remain; XML dimensions are metadata, not recomputed cell counts.

## Supplementary table S3_Method comparison in cell type annotation.xlsx

Original path: `paper/submission_v16/Supplementary table S3_Method comparison in cell type annotation.xlsx`

SHA-256: `b6140dfa828777ed51e9d2330715e4754d1b96a38cca2707c752b60884e0fe6e`

Bytes: 12374640.

Container CRC: passed. Workbook sheets: `Supplementary table S3_Method c`.

Declared XML sheet dimensions: `xl/worksheets/sheet1.xml: A1:Z92410`. These are metadata, not newly computed cell counts.

Remaining: Full spreadsheet reader, endpoint/content validation and submission cross-reference check remain; XML dimensions are metadata, not recomputed cell counts.

## Supplementary table S4_PTC_cell type composition.xlsx

Original path: `paper/submission_v16/Supplementary table S4_PTC_cell type composition.xlsx`

SHA-256: `d5277be9031982248529e620532f7f523120b49844ebe70570519708f8742622`

Bytes: 11588.

Container CRC: passed. Workbook sheets: `Supplementary Table 4`.

Declared XML sheet dimensions: `xl/worksheets/sheet1.xml: A1:F27`. These are metadata, not newly computed cell counts.

Remaining: Full spreadsheet reader, endpoint/content validation and submission cross-reference check remain; XML dimensions are metadata, not recomputed cell counts.

## Content-specific acceptance

| Original item | Required revision checks |
|---|---|
| Figure S1 | Render the original PDF/PNG and revised submission figure; verify legible labels, method names, preprocessing/input provenance and captions. A separated UMAP picture does not establish superior accuracy or preserved original-space density. Bind any empirical claim to actual partition/terminal metrics. |
| Table S1 | Check gene lists, tissue/context and assay provenance; distinguish AllTissues seven-tissue union from all-human/all-tissue catalogues. The first workbook sheet is named `Supplementary Table 2` despite the S1 filename; correct only in a versioned revision copy and document the change. Map this table to the exact marker configurations used in current runs. |
| Table S2 | Preserve the original final annotations; verify sample/barcode identity, headers and endpoint. Do not rename TCR support as independent T-cell-subtype validation. Reconcile its historical labels independently of S3. |
| Table S3 | Preserve native predictions, original flags and historical SignacX zero-T row. Separate the paper Table2 metric reconstruction, literal S3 evaluation and new strict-T comparisons. Verify all-cell denominators, ontology mappings and historical Accuracy disclosure. |
| Table S4 | Trace each composition value to its S2/S3 endpoint, full denominator and sample/patient identity; retain paired patient structure and contextual limits. |

Final delivery needs a supplement index, direct readable previews where useful,
machine-readable tables, consistent numbering and tested local/packaged links.
The current inventory does not claim this final submission-level acceptance.
