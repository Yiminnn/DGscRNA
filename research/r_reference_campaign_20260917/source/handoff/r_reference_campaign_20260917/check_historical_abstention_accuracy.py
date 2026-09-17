"""Test a specific metric-definition hypothesis against the reported DG Accuracy entry."""
import os,json
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import pandas as pd
    base=ROOT/'results/hvg_ptc_20260916_v1/ptc_paper_baseline'
    ref=pd.read_csv(base/'paper_baseline_reference.csv.gz',keep_default_na=False).set_index('cell_id')
    pairs=pd.read_csv(base/'original_DG_binary_pairs_for_R.csv.gz').set_index('cell_id').loc[ref.index]
    target=pd.read_csv(base/'paper_Table2_target_comparison.csv')
    target=target[target.method.eq('DG')&target.metric.eq('Accuracy')]
    known=~ref.paper_final_native.isin(['Unknown','Undecided','No_Annotation'])
    correct=pairs.truth.eq(pairs.prediction)
    rows=[]
    for row in target.itertuples(index=False):
        mask=pairs.index.notna() if row.scope=='Overall' else pairs.scope.eq(row.scope)
        accuracy_error=float((correct[mask]&known[mask]).mean())
        accuracy_called=float(correct[mask&known].mean())
        rows.append(dict(scope=row.scope,paper_accuracy=row.paper,n_cells=int(mask.sum()),
          historical_binary_mapping_accuracy=float(correct[mask].mean()),unknown_fraction=float((~known[mask]).mean()),
          unknown_as_error_accuracy=accuracy_error,called_only_accuracy=accuracy_called,
          unknown_as_error_matches_4dp=round(accuracy_error,4)==row.paper,
          called_only_matches_4dp=round(accuracy_called,4)==row.paper))
    report=dict(job=os.environ['SLURM_JOB_ID'],hypothesis='Manuscript Accuracy used an alternative Unknown policy while F1 used the historical binary mapping.',
      rows=rows,all_unknown_as_error_rows_match=all(r['unknown_as_error_matches_4dp'] for r in rows),
      all_called_only_rows_match=all(r['called_only_matches_4dp'] for r in rows),
      limitation='Numerical agreement alone would not recover the original computation script; no archived labels are changed.')
    (OUT/'verification/historical_accuracy_unknown_policy_hypothesis.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':run()
