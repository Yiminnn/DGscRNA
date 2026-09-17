"""Compare terminal original-route refits with archived cells and Table 2."""
from pathlib import Path
import os,sys,json,hashlib
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
BASE=ROOT/'results/hvg_ptc_20260916_v1'
OUT=BASE/'ptc_paper_baseline'

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import numpy as np
    import pandas as pd
    from sklearn.metrics import f1_score,roc_auc_score,accuracy_score
    flavor=sys.argv[1]
    assert flavor in ['full','full_parallel','marker_union']
    parent=OUT/({'full':'replay_selected_routes','full_parallel':'replay_selected_routes_full_parallel',
                 'marker_union':'replay_selected_routes_marker_union'}[flavor])
    dest=OUT/('evaluation_'+flavor);dest.mkdir(exist_ok=True)
    ref=pd.read_csv(OUT/'paper_baseline_reference.csv.gz',keep_default_na=False).set_index('cell_id')
    y=pd.read_csv(OUT/'original_DG_binary_pairs_for_R.csv.gz',keep_default_na=False).set_index('cell_id').loc[ref.index]
    T_names=set(json.loads((OUT/'historical_T_names_from_vignette.json').read_text()))
    def general(x):
        if x in {'Unknown','Undecided','No_Annotation'}:return x
        fields=x.split('+')
        if fields[0]=='NCOMMREFF':return fields[1]
        if fields[0]=='cancer':return '+'.join(fields[3:])
        if fields[0] in ['CellMarker_normal','CellMarker_cancer']:return '+'.join(fields[3:])
        return '+'.join(fields[2:]) if len(fields)>1 else x
    frames=[];manifest=[];parities=[]
    for route,group in [('NMT_Thyroid_Seurat_none','NMT'),('TTU_Pubmed_UMAPHDBSCAN_mean','TTU')]:
        d=parent/route
        assert (d/'TERMINAL_COMPLETE').exists(),str(d)
        m=json.loads((d/'terminal_DL/training_manifest.json').read_text())
        f=pd.read_csv(d/'terminal_predictions.csv.gz',keep_default_na=False).set_index('cell_id')
        assert f.index.is_unique and set(f.index)==set(ref.index)
        f=f.loc[ref.index]
        selected=ref.group.eq(group)
        f=f.loc[selected].copy()
        native=ref.loc[selected,'paper_final_native']
        f['historical_native']=native
        f['route']=route
        f['initial_known']=~f.initial_native.isin(['Unknown','Undecided','No_Annotation'])
        for stage in ['initial_native','terminal_native_090','terminal_native_070_diagnostic']:
            f[stage+'_exact']=f[stage].eq(native)
            f[stage+'_general']=f[stage].map(general)
            f[stage+'_historical_broad_T']=f[stage+'_general'].isin(T_names).astype(int)
        known_mismatch=f.initial_known & f.initial_native.ne(native)
        manifest.append(dict(route=route,group=group,n_cells=len(f),dl_status=m['dl_status'],
            n_initially_known=int(f.initial_known.sum()),n_initial_known_conflict_with_archive=int(known_mismatch.sum()),
            n_terminal_exact=int(f.terminal_native_090_exact.sum()),
            n_terminal_mismatch=int((~f.terminal_native_090_exact).sum()),
            n_terminal_unresolved=int(f.terminal_native_090.isin(['Unknown','Undecided','No_Annotation']).sum())))
        f.loc[~f.terminal_native_090_exact].to_csv(dest/(group+'_terminal_mismatches.csv.gz'))
        frames.append(f)
        if flavor in ['full','full_parallel']:
            other=OUT/'replay_selected_routes_marker_union'/route
            assert (other/'TERMINAL_COMPLETE').exists()
            fast=pd.read_csv(other/'terminal_predictions.csv.gz',keep_default_na=False).set_index('cell_id').loc[ref.index]
            full=pd.read_csv(d/'terminal_predictions.csv.gz',keep_default_na=False).set_index('cell_id').loc[ref.index]
            a=pd.read_csv(other/'DEG_marker_union.csv.gz').set_index(['cluster','gene']).sort_index()
            b=pd.read_csv(d/'DEG_full_integrated2000.csv.gz').set_index(['cluster','gene']).sort_index()
            assert a.index.is_unique and b.index.is_unique and a.index.isin(b.index).all()
            b=b.loc[a.index]
            cols=['p_val','avg_log2FC','pct.1','pct.2','p_val_adj']
            delta=float(np.max(np.abs(a[cols].to_numpy()-b[cols].to_numpy())))
            parity=dict(route=route,n_tested_rows=len(a),max_absolute_statistic_delta=delta,
                initial_exact=int(full.initial_native.eq(fast.initial_native).sum()),
                terminal_exact=int(full.terminal_native_090.eq(fast.terminal_native_090).sum()),n_cells=len(full))
            parities.append(parity)
            assert delta<=1e-12 and parity['initial_exact']==parity['terminal_exact']==len(full)
    final=pd.concat(frames).loc[ref.index]
    final['paper_tissue']=y.scope
    final['TCR_original_any_contig']=y.truth
    final.to_csv(dest/'combined_selected_terminal_predictions.csv.gz')
    rows=[]
    scopes=[('Overall',final.index)]+list(final.groupby('paper_tissue').groups.items())+list(final.groupby('group').groups.items())
    for stage in ['initial_native','terminal_native_090','terminal_native_070_diagnostic']:
        for scope,ids in scopes:
            yy=final.loc[ids,'TCR_original_any_contig'].astype(int)
            pred=final.loc[ids,stage+'_historical_broad_T'].astype(int)
            rows.append(dict(stage=stage,scope=scope,n=len(ids),
                F1_source_default_positive0=f1_score(yy,pred,pos_label=0),
                F1_T_positive1=f1_score(yy,pred,pos_label=1),AUC_from_binary_calls=roc_auc_score(yy,pred),
                accuracy=accuracy_score(yy,pred)))
    metrics=pd.DataFrame(rows);metrics.to_csv(dest/'selected_route_metrics.csv',index=False)
    targets=pd.read_csv(OUT/'paper_Table2_historical_source_comparison.csv')
    targets=targets[targets.method.eq('DG')].copy()
    cols={'F1 score':'F1_source_default_positive0','AUC-ROC':'AUC_from_binary_calls','Accuracy':'accuracy'}
    targets['rerun_terminal090']=[metrics.loc[metrics.stage.eq('terminal_native_090') & metrics.scope.eq(r.scope),cols[r.metric]].item() for r in targets.itertuples()]
    targets['rerun_minus_archived']=targets.rerun_terminal090-targets.reconstructed
    targets['rerun_matches_paper_4dp']=[round(r.rerun_terminal090,4)==r.paper for r in targets.itertuples()]
    targets.to_csv(dest/'paper_Table2_terminal_rerun_comparison.csv',index=False)
    pd.DataFrame(manifest).to_csv(dest/'terminal_retraining_agreement.csv',index=False)
    if parities:pd.DataFrame(parities).to_csv(dest/'marker_union_full_gene_parity.csv',index=False)
    summary=dict(job=os.environ['SLURM_JOB_ID'],flavor=flavor,routes=manifest,
        terminal_exact_cells=sum(r['n_terminal_exact'] for r in manifest),n_cells=len(final),
        DG_table_rows_matching_paper=int(targets.rerun_matches_paper_4dp.sum()),
        original_accuracy_source_unresolved=True,
        status='reconstruction_measured; do not automatically mark full paper gate passed',
        full_gene_marker_union_parity=parities)
    (dest/'manifest.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(pd.DataFrame(manifest).to_string(index=False),flush=True)
    print(targets[['scope','metric','paper','reconstructed','rerun_terminal090','rerun_minus_archived','rerun_matches_paper_4dp']].to_string(index=False),flush=True)
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':run()
