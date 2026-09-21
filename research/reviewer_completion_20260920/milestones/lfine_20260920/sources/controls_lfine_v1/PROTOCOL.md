# B: Lfine compatible-target-set metrics

Frozen post-fit evaluation of the completed three-pilot GBM neighbor and learning
controls. No fitting, label changes, checkpoint selection, or new optimum search.
The original L1 evaluation remains archived separately and contributes no numerical
performance values to this output.

## Endpoint and mapping

Use the exact current compact evaluator, v5 helper implementation and frozen native
panel semantic map, with source and truth SHA-256 pins in `freeze.json`. The copied
compact provider is imported only for `initialize`; none of its output-writing
entrypoints is called. Its v5 evaluator and AST-extracted original `grid_sets.py`
helpers determine label composition and metrics. Native panel predictions remain
unchanged. Semantic parents map to compatible Lfine target sets through frozen
`LFINE_PREFIX`; this is not one-to-one fine-type prediction.

For each sample, score classes with at least 20 reference cells, excluding Other
and nan. Build compatible target sets from every observed reference class and keep
all cells in TP/FP/FN. Unknown, Undecided and unmappable calls get no compatible
credit. Do not infer a cell's predicted subtype from its true subtype. The source
map's internal L1 column is ontology provenance only. No L1 metric is output.

CM2_glioma_other retains the exact historical v5 mapping. The prespecified
CM2_primary_all_context map extends the original marker vocabulary; it was not
learned from the new outcome data. All panel mappings and target prefixes are
exported for inspection. Missing/invalid terminal states cannot become marker-only
results. Every prediction file, terminal provenance, truth order and known-label
retention is checked; source hashes are checked again after evaluation.

## Fixed scope

- Samples: TKU4163, NL022, SN040; budgets: HVG2000 and HVG5000.
- 66 unique neighbor configurations, two fixed markers, two thresholds: 264 rows.
- 18 existing learning fits, epochs 5/10/20/30, two thresholds: 144 rows.
- 408 total terminal-only metric rows. Threshold 0.90 is primary; 0.70 is sensitivity.
- Markers: CM2_glioma_other and CM2_primary_all_context, mean seed cutoff. Learning
  uses the primary marker only; backup markers are not mixed into that series.
- All epoch checkpoints and seeds 0/1/42 remain reported, with no outcome-based
  selection. Training validation curves measure held-out marker pseudo-labels.

## Pilot gate

Only 24 terminal090 rows are authorized initially: 18 default-neighbor endpoints
(three samples, two budgets, PCA30_SNN/UMAP2_SNN/UMAP2_HDBSCAN_R at k20/n30) and six
learning seed42/epoch10 checkpoints. Compare all 12 compact SCORE_FIELDS with the
saved compact primary-marker results, tolerance 1e-12 and matching NA states.
Require exact native initial/final prediction and cell-order identity against
the compact original-R terminal artifacts as well. Independently recompute each
scored class's compatible TP/FP/FN and macro-F1.

Full evaluation requires a recorded parent authorization after the pilot passes.
The pilot does not establish sensitivity conclusions or uniform optimality.

## Compact notebook integration

Do not modify the canonical notebook from this package. Supply a later, concise
HVG/parameter sensitivity result: at most one Lfine neighbor figure and a compact
checkpoint/coverage table, retaining all raw rows and per-class counts externally.
These HVG-only controls do not supply the all-gene main figure. Do not reuse old
L1 plots or conclusions under an Lfine heading. Unfavorable comparisons remain
in the complete auditable results and are considered when assessing any claim.
