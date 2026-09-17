"""Trace disagreements already observed; no fitting, tuning, or label modification."""
from pathlib import Path
import os,json
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
BASE=ROOT/'results/hvg_ptc_20260916_v1'
OUT=BASE/'ptc_paper_baseline'

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import numpy as np,pandas as pd
    d=pd.read_csv(OUT/'evaluation_marker_union/combined_selected_terminal_predictions.csv.gz',keep_default_na=False)
    rows=[];transitions=[]
    unknown=['Unknown','Undecided','No_Annotation']
    for group,g in d.groupby('group'):
        for stage in ['terminal_native_090','terminal_native_070_diagnostic']:
            old=g.historical_native;new=g[stage]
            bad=new.ne(old)
            assert (~g.initial_known[bad]).all()
            broadnew=g[stage+'_historical_broad_T']
            hist=g.historical_native.map(lambda x:x.split('+',1)[1] if x.startswith('NCOMMREFF+') else '+'.join(x.split('+')[3:]) if x.startswith('cancer+') else x)
            histT=hist.isin(json.loads((OUT/'historical_T_names_from_vignette.json').read_text())).astype(int)
            rows.append(dict(group=group,stage=stage,n_cells=len(g),n_mismatch=int(bad.sum()),
                old_known_new_unknown=int((bad & ~old.isin(unknown) & new.isin(unknown)).sum()),
                old_unknown_new_known=int((bad & old.isin(unknown) & ~new.isin(unknown)).sum()),
                changed_known_cell_type=int((bad & ~old.isin(unknown) & ~new.isin(unknown)).sum()),
                changed_binary_T_NK_status=int(broadnew.ne(histT).sum()),
                new_unknown=int(new.isin(unknown).sum()),old_unknown=int(old.isin(unknown).sum())))
            for (a,b),count in g.loc[bad].groupby(['historical_native',stage]).size().items():
                transitions.append(dict(group=group,stage=stage,historical_native=a,new_native=b,n_cells=int(count)))
    pd.DataFrame(rows).to_csv(OUT/'DL_discrepancy_categories.csv',index=False)
    pd.DataFrame(transitions).to_csv(OUT/'DL_discrepancy_transitions.csv',index=False)
    # Full-cohort original Pubmed branch exists independently in the extended CSV.
    archive=BASE/'ptc_recovery/archive/tcr/rawdata'
    e=pd.read_csv(archive/'data_with_validation+3_cell_types 2.csv',keep_default_na=False)
    meta=pd.read_csv(archive/'metadata.txt',sep='\t')
    samplemap=dict(zip(meta.sc_ID,meta.Sample))
    e['cell_id']=e['sample.name'].map(samplemap)+'_'+e.X.str.extract(r'([ACGT]{12,})',expand=False)
    e=e.set_index('cell_id')
    branch='hdbscan.UMAP_clusters_NCOMMREFF_mean_DGCyTOF_General'
    fresh=pd.read_csv(OUT/'replay_selected_routes_marker_union/TTU_Pubmed_UMAPHDBSCAN_mean/terminal_predictions.csv.gz',keep_default_na=False).set_index('cell_id')
    e=e.loc[fresh.index]
    known=fresh.initial_native.ne('Undecided')
    initial=fresh.initial_native.str.replace(r'^NCOMMREFF\+','',regex=True)
    conflicts=known & initial.ne(e[branch])
    report=dict(job=os.environ['SLURM_JOB_ID'],Pubmed_full_cohort_n_cells=len(fresh),
        Pubmed_full_cohort_n_initial_known=int(known.sum()),
        Pubmed_full_cohort_initial_known_conflicts=int(conflicts.sum()),
        no_new_models_fitted=True,no_reference_labels_changed=True,
        source_DL_threshold=0.9,manuscript_threshold_diagnostic_only=0.7,
        historical_model_initialization='not recorded in recovered source or checkpoint; cause of exact terminal differences cannot be uniquely assigned')
    (OUT/'DL_discrepancy_trace_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    print(pd.DataFrame(rows).to_string(index=False),flush=True)
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':run()
