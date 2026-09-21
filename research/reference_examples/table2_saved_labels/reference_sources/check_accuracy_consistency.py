"""Check the recovered single-class F1 against the paper's Accuracy row.

For the same binary calls, F1=2TP/(2TP+E), Accuracy=(TP+TN)/N,
and TP<=TP+TN imply F1<=2*Accuracy/(1+Accuracy).
This algebra cannot identify the unknown source of the Accuracy row.
"""
import os,json
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/ptc_paper_baseline'

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import pandas as pd
    table=pd.read_csv(OUT/'paper_Table2_historical_source_comparison.csv')
    overall=table[table.scope.eq('Overall')]
    rows=[]
    for method,g in overall.groupby('method'):
        f=g.loc[g.metric.eq('F1 score'),'paper'].item()
        a=g.loc[g.metric.eq('Accuracy'),'paper'].item()
        lower_f=f-0.00005;upper_a=a+0.00005
        min_accuracy=lower_f/(2-lower_f)
        rows.append(dict(method=method,paper_F1=f,paper_accuracy=a,
            minimum_accuracy_allowed_by_F1_with_rounding=min_accuracy,
            paper_accuracy_upper_rounding_bound=upper_a,
            compatible_same_binary_prediction_vector=bool(upper_a>=min_accuracy),
            recalculated_accuracy=g.loc[g.metric.eq('Accuracy'),'reconstructed'].item()))
    out=pd.DataFrame(rows);out.to_csv(OUT/'paper_accuracy_same_endpoint_consistency.csv',index=False)
    report=dict(job=os.environ['SLURM_JOB_ID'],n_methods=len(out),
        n_incompatible=int((~out.compatible_same_binary_prediction_vector).sum()),
        applicable_scope='Source-recovered one-class binary F1 and ordinary accuracy on the same predictions and denominator',
        bound='F1 <= 2*Accuracy/(1+Accuracy), equivalently Accuracy >= F1/(2-F1)',
        conclusion='Table 2 Accuracy cannot be the ordinary accuracy of the same binary evaluation. Its original source or intended different endpoint is not recovered; do not alter predictions to target it.')
    (OUT/'paper_accuracy_consistency_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    print(out.to_string(index=False));print(json.dumps(report,indent=2))

if __name__=='__main__':run()
