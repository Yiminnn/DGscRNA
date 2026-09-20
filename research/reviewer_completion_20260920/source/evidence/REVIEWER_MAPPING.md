# GBM evidence supplement: reviewer mapping

Scope is the already frozen GBM cohort and saved terminal predictions. This work
does not change the existing annotation endpoint, select a new marker library,
or add model fits. Execution source is `gbm_evidence.py`, SLURM entrypoint is
`run.sbatch`, and accounting verification is `verify_evidence.py`.

| Requirement | New evidence | Completion boundary |
|---|---|---|
| R3-S1: cross-tool consistency and clustering metrics | Six actual tools; pairwise agreement, ARI/NMI/FMI; marker-only four-tool and all-six consensus; accuracy versus author L1 in agreement/disagreement strata; full cells, called cells and supported-label denominators; patient aggregation | GBM component. Annotation-induced partitions are clearly distinguished from upstream clustering. Author labels remain a reference annotation, not independent biological truth. |
| R3-S6: expression of selected markers | NL022 true normalized-RNA dotplot/violin; all selected genes and original sources; full-HVG R normalization parity; both fixed and training-patient-selected contexts on unchanged coordinates | GBM component. Selected library is frozen from the existing patient-label-heldout selection. Expression is consistency evidence, not independent validation. PTC counterpart remains in its later work package. |
| R2-8: marker performance / interpretation | Concrete selected genes, assays and label endpoints; existing cohort marker-performance distributions retained | Complements existing performance distributions; does not replace or re-rank them. |
| PI recent single-sample display request | Frozen NL022 author labels, fixed-marker initial and terminal calls, selected-marker initial and terminal calls, Unknown on identical native-R UMAP coordinates | Expression is not used to revise the author or terminal labels. |

Outputs are placed under
`results/hvg_ptc_20260916_v1/reviewer_completion_20260920/evidence/`.
The webpage should consume that directory's `status.json`. Only a completed
manifest, successful denominator verification and figure inspection justify
marking this bounded GBM subtask complete. Resource-ledger work is a separate
part of F and is not covered by this supplement. Historical SignacX 0 T cells
and every PTC endpoint are untouched.
