# Source-grounded Methods replacement draft

Status: proposed replacement text, 2026-09-20. The original manuscript is unchanged. Bracketed fields are unresolved or depend on ongoing experiments. Source positions `Pxxxx` refer to the extracted original V16 DOCX paragraphs in `original_v16.paragraphs.txt`, not submission page numbers. These source corrections do not establish that every reference setting is optimal.

## Framework and input data

DG-scRNA combines upstream cell partitioning, context-dependent marker annotation and expression-based multilayer-perceptron refinement. The classifier consumes expression features; graph edges are not classifier inputs, and the implementation does not perform message passing or joint graph-neural-network learning. The workflow is discriminative. Its original hybrid implementation uses R for preprocessing, integration, clustering and marker scoring, and Python for terminal refinement.

For the PTC historical benchmark, the delivered GEX inputs contain 110,497 sample-specific barcodes. The preserved downstream annotation endpoint contains 92,404 cells. These counts refer to different stages. [INSERT the verified barcode-level attrition table, including explicitly unavailable intermediate reasons.] The original S2 and S3 annotations represent distinct final endpoints and are retained separately. Original delivered files, reconstruction outputs and new refits are identified by source hashes and run identifiers.

Sources: V16 P0037-P0038, P0315, P0563, P0610-P0615, P0657; `handoff/original_rmd_workflow_audit_20260920/rmd_source_audit.md`; `results/hvg_ptc_20260916_v1/geo_submission_v1/PTC_GEO_submission_20260916/README.md`; `results/hvg_ptc_20260916_v1/ptc_paper_baseline/home_original_sup_manifest.json`.

## Preprocessing and feature roles

The recovered R helper retains cells with `nFeature_RNA > 200` and mitochondrial percentage `< 15`, calculated from genes with the `^MT-` prefix. Seurat LogNormalize uses a scale factor of 10,000. Variable features are selected by the VST method with `nfeatures=2000`. Before doublet detection, the helper scales the data while regressing mitochondrial percentage and computes 30 principal components. DoubletFinder uses a pK sweep, selects the maximal BCmetric, and applies PCs 1-10, pN=0.25 and an expected doublet count of round(0.075 × the current number of cells), retaining Singlet calls. [INSERT independently verified historical package version; do not equate the recovery environment with the original run.]

The helper documents these thresholds, but available intermediate historical cell membership must still be checked before assigning every excluded barcode to a particular filter. Gene filtering at object creation is distinct from cell filtering and downstream feature selection. No upper detected-gene cutoff is added to the reference description without source evidence.

For multi-sample CCA integration, `SelectIntegrationFeatures`, `FindIntegrationAnchors` and `IntegrateData` are applied. The saved historical eight-sample PTC object records 2,000 integration features, LogNormalize/CCA and dimensions 1-30. The integrated assay is scaled and PCA30 is computed; the integration-stage scaling call does not repeat the explicit mitochondrial-regression argument used in sample-level preprocessing. The archived object used eight-sample joint CCA; new analyses specified by the current protocol integrate NMT (MT-1/2 and N-1/2) and TTU (TU-1/2 and T-1/2) separately. These are distinct named analyses.

The manuscript will distinguish (i) sample-level variable genes, (ii) integration features, (iii) representation genes, (iv) the assay and eligible genes used for differential expression and marker scoring, and (v) the terminal classifier's features. In the single-sample GBM reference, marker scoring uses eligible RNA genes while representation and DL use their declared selected features. In the historical integrated PTC branch, the saved integrated matrix contains 2,000 features and the active assay determines the DEG/scoring universe. The ScaleData parameter `min.cells.to.block=3000` is a computational blocking parameter, not an HVG count. The recovered reference does not support the manuscript's TMM/edgeR normalization stage.

Sources: archived `results/hvg_ptc_20260916_v1/ptc_recovery/archive/tcr/ptc_val/scripts/DGscRNA-Share/R/source.R:94-151,204-226,338-345`; `results/hvg_ptc_20260916_v1/ptc_recovery/inventory_workspace/object_01.commands.txt`; `objects.json` in the same folder; V16 P0599-P0601, P0664.

## Representation and partitioning

The reference UMAP branch applies uwot UMAP to PCs 1-30, with two output dimensions, 30 neighbors, cosine distance, minimum distance 0.3 and the saved representation seed 42. The recovered framework also contains PCA- and UMAP-based SNN and HDBSCAN candidates. The reference SNN clustering resolution is 0.5; R HDBSCAN uses minPts=50. The historical PTC selected contexts are NMT with CellMarker Thyroid, PCA-SNN and no density cutoff, and TTU with Pubmed34663816-derived markers, UMAP-HDBSCAN and the mean density cutoff. Both then use terminal refinement. Thus, UMAP-HDBSCAN is not the historical selected route for both groups.

The new seven-representation comparison is a separately specified workflow experiment. It uses the same scaled-HVG matrix as the representation input, whereas the original UMAP anchor uses PCA30 first. These distinct inputs and paths will remain labeled. [INSERT completed, validated comparison coverage and parameters.] Partition metrics and terminal annotation metrics will not be combined under an unlabeled accuracy heading.

Sources: archived R helper lines 235-330; saved Seurat command records; `results/hvg_ptc_20260916_v1/ptc_recovery/inventory/notebook_sources/06_ptc_paper.txt:208-239`; approved `handoff/reviewer_completion_plan_20260920/PLAN.md`, work package A.

## Marker density and initial labels

For cluster c and candidate cell type t, let M_t be the full marker panel and D_c the cluster DEG list. The implemented score is the sum of cluster average log2 fold changes for genes in the intersection of M_t and D_c with average log2 fold change greater than 1, divided by the full panel length |M_t|. A marker panel of length at most one receives the implemented 0.8 multiplier. The highest-scoring cell type is assigned to the cluster; a tie produces Undecided. The named cutoff conditions are none, 0.5 and mean. The mean condition compares a cluster's maximum score against the mean of the maximum scores across clusters. The numerical cutoff can therefore change when the partition changes.

This is a cluster-DEG score, not a mean-expression or alpha-weighted Jaccard mixture. The source does not implement the manuscript's alpha parameter or epsilon-based iterative convergence rule. The new no-clustering control necessarily replaces this seed-construction rule and will describe that replacement explicitly.

Sources: archived R helper lines 378-430; `results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/GBM_full_summary/METHODS_CORRECTIONS.md`; V16 density/formula paragraphs preceding P0657. Formula objects not represented in the plain-text extraction must also be checked in the final DOCX/PDF.

## Terminal DL/refinement

Initially assigned cells provide the training labels, and only Undecided cells are eligible for new assignments. The preserved network has hidden layers of 256 and 128 units with LeakyReLU activations and a Softmax output. The reproduction retains the original Softmax-output/CrossEntropyLoss combination, Adamax at learning rate 0.001, a batch size of 256, ten epochs and the original non-shuffled training order. The archived 90/10 split uses random state 42. A global initialization seed is not established by that split setting; new diagnostic runs record initialization seeds separately.

The class vocabulary follows the original implementation, including its representation of Undecided when present. Initially known labels are retained. Predictions for unresolved cells are accepted when the original rounded probability rule, round(probability, 4) ≥ 0.90, is met; lower-confidence cells remain Unknown. The 0.70 threshold is a separately reported sensitivity condition. When all cells are already assigned, refinement is a no-op; when no usable labeled cells exist, a learned-refinement result is unavailable. These states are reported rather than counted as successful fitted models.

The new training diagnostic records held-out pseudo-label validation at every epoch and prespecified checkpoints through 30 epochs. It does not establish independent biological validation, and no test truth is used for early stopping. [INSERT the validated curves, checkpoints, scopes and any legal no-training states.]

Sources: archived `R/source.py:47-69,119-147,165-209,386-435`; V16 P0657, P0666; `handoff/reviewer_completion_plan_20260920/PLAN.md`, work package B. The original weights remain unavailable; refit equivalence is not represented as bitwise recovery of an unavailable historical training run.

## Evaluation, model selection and historical metrics

DG-scRNA annotation evaluation uses the terminal DL/refinement output. Marker-only outputs, partition evaluations and no-op or invalid training states are identified separately. New multiclass evaluation uses a frozen label vocabulary, per-class metrics, macro-F1, weighted-F1 and coverage, with Unknown and unsupported classes retained under the declared all-cell denominator rule. Patient-level splitting and paired inference keep all samples from the same patient together. Marker and parameter selection uses training patients only. Shared unlabeled integration is named transductive rather than fully inductive generalization. [INSERT cohort-specific fold, class and selection manifests.]

The original PTC historical F1 reconstruction uses the recovered original reference definition and MLmetrics' class-0-positive convention; historical AUC from binary calls is distinguished from probability-based AUC. Strict T-positive, literal S3 and alternate TCR-definition results are separate named evaluations. The historical Accuracy definition has not been recovered and will not be presented as a reproduced ordinary-accuracy value. TCR positivity supports T-cell detection but cannot validate every non-T class or T-cell subtype. Its absence is not independent proof of non-T identity.

The historical SignacX table remains None with zero predicted T cells, as originally reported. The populated native non-T predictions, the unresolved reason for zero T calls, and later nonzero predictions are reported in separate provenance and revision-comparison tables. scType retains its actual historical results. Headline fair comparisons require matched information and input conditions; preserving a historical table does not establish those conditions.

Sources: V16 P0128, P0662, P0666-P0671; `handoff/reviewer_completion_plan_20260920/SCTYPE_AUDIT.md`; `results/hvg_ptc_20260916_v1/ptc_paper_baseline/MLmetrics_actual_definitions.txt`; `paper_Table2_historical_source_comparison.csv` and `paper_accuracy_same_endpoint_consistency.csv` in the same folder.

## Text that must remain unresolved until original records are available

Do not silently choose >80% or >85% viability, 70-µm or 40-µm strainer, a definitive historical DoubletFinder version, or missing library chemistry. Replace those entries only after an original protocol/run record is attached to SOURCE_GAPS.md. For the six available original GEX reports, the chemistry is recorded as Single Cell 5′ PE; the original manuscript's 3′ assertion requires correction or documented reconciliation. The available original VDJ reports say Cell Ranger 3.0.1, whereas V16 P0662 says 7.1. Do not use the current software environment as retrospective evidence of either value.
