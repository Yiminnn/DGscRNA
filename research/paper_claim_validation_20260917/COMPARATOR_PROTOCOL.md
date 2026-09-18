# Frozen comparator interpretation

All fitting, audits, evaluation and plots run through SLURM. This protocol records the shared information and valid scope before inspecting comparator accuracy rankings. It does not claim a completed benchmark.

## Primary matched-marker comparison

- Compare terminal DG-scRNA, official-audited scType, and official scCATCH at the same original native-R HVG2000 -> PCA30 -> UMAP2 -> R-HDBSCAN50 partition. Include every cell, including noise/Unknown/unmapped calls.
- Give each method all16frozen marker libraries. Primary retrospective patient-heldout marker selection excludes CARE_TME, BrainAtlas112 and UNION_all because they contributed to author-label construction. Complete concordance results retain them.
- Each cluster-based method has3predeclared method-specific thresholds: DG density none/mean/0.5, scType cluster score/cell count0.125/0.25/0.5, scCATCH Wilcoxon p0.01/0.05/0.1. These are comparable selection opportunities, not numerically identical parameters.
- scType also runs all24native-R clustering configurations because its exact linear aggregation is inexpensive. The DG24-route/feature selection and scType24-route/feature selection are reported separately from the matched-partition primary comparison. A tuned24-configuration DG result must not be compared to fixed-partition scCATCH as if tuning opportunities were identical.
- scCATCH uses its original pairwise-cluster DEG method1, logFC0.25 and expression proportion0.25. Gene-level tests can be reused across frozen marker libraries only after exact native-function parity. Its narrower geometry scope is a computational budget decision recorded without using prediction accuracy. All16libraries and3cutoffs remain in scope.
- SCINA is a per-cell model and does not consume a clustering partition. It receives the same marker libraries and full normalized RNA, original defaults and overlap-removal0/1. Default overlap-removal1 remains its anchor. Constant/empty signatures, numerical errors and any audited numerical boundary guard are reported explicitly. Numerical errors are not valid zero-performance arms.
- Marker native labels are preserved until final mapping. scCATCH tied labels may collapse to one L1 parent only when all exact native tied panels map to that same parent; otherwise they remain unmappable. Test labels do not resolve ties.

## Reference-information comparison

- SingleR2.8.0 uses the same five patient folds; all samples of each heldout patient are excluded from reference construction. Training labels come only from primary97samples. Reference cells are selected by a frozen barcode hash, at most100per training-patient x L1 class; rare-class cells are retained.
- SingleR uses its built-in reference aggregation, exact neighbors and fine tuning. Default pruned per-cell calls are primary; unpruned calls are a sensitivity. Gene availability can be read from heldout input metadata, but heldout labels are loaded only after prediction.
- This is a supervised training-patient reference condition and is shown separately from marker-only information. It must not be described as needing the same prior information as DG-scRNA.
- Pretrained scDeepSort, if its published CPU environment is recovered/installed successfully, is a separate atlas-trained GNN condition. Unsupported malignant/other classes remain in the denominator. Dependency errors are not biological failures.

## Evaluation

Use frozen L1 mapping, present-class macro-F1 with every cell in the denominator, fixed11secondary, per-class confusion, Unknown/coverage, and patient-weighted means. Select libraries/thresholds on training patients only. Report paired heldout-patient differences with patient bootstrap intervals and multiplicity correction. These are retrospective folds in an explored cohort; no never-seen-cohort claim is made.

Official references: [scType](https://github.com/IanevskiAleksandr/sc-type), [scCATCH](https://github.com/ZJUFanLab/scCATCH), [SCINA1.2.0](https://search.r-project.org/CRAN/refmans/SCINA/html/SCINA.html), [SingleR advanced workflow](https://bioconductor.org/books/release/SingleRBook/advanced-options.html).
