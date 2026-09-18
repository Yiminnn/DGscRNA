"""QC and refinement-error associations; no doublet or novel-cell-type inference."""
import os
from common import OUT, require_slurm, checked, complete, sha, write_json, utc

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from prediction_helpers import read_native,map_labels,UNKNOWN

    assert checked(OUT/'comparison_summary')
    dest=OUT/'unknown_summary';dest.mkdir(exist_ok=True)
    cohort=pd.read_csv(OUT/'protocol/cohort.csv')
    cv=pd.read_csv(OUT/'comparison_summary/patient_heldout_results.csv',dtype={'cutoff':str})
    cv=cv[(cv.method=='DG-scRNA')&(cv.cohort=='primary97')]
    rows=[];errors=[];sources={}
    for s in cohort[cohort.primary].itertuples():
        truth=pd.read_csv(OUT/'evaluation_inputs'/s.sample/'truth.csv.gz',dtype=str,keep_default_na=False)
        qc=pd.read_csv(OUT/'evaluation_inputs'/s.sample/'qc.csv.gz')
        assert list(qc.CellID)==list(truth.cell_id)
        chosen=cv[cv.patient==s.patient].iloc[0]
        configurations=[('fixed_glioma','CM2_glioma_other','mean'),('training_patient_selected',chosen.library,chosen.cutoff)]
        for condition,library,cutoff in configurations:
            p,source=read_native('DG-scRNA',s.sample,'hvg2000','UMAP2_HDBSCAN_R',library,cutoff)
            assert list(p.cell_id)==list(truth.cell_id)
            sources[str(source.relative_to(OUT))]=sha(source)
            mapped=map_labels('DG-scRNA',library,p.prediction)
            seed=map_labels('DG-scRNA',library,p.initial)
            undecided=p.initial.eq('Undecided').to_numpy()
            unknown=p.prediction.isin(UNKNOWN).to_numpy()
            category=np.where(unknown,'Unknown',np.where(~undecided,'retained_seed',np.where(mapped==truth.L1,'new_correct','new_incorrect')))
            data=qc.copy();data['truth']=truth.L1;data['category']=category;data['unknown']=unknown
            for group,indices in data.groupby(['truth','category']).groups.items():
                take=data.loc[indices]
                rows.append(dict(sample=s.sample,patient=s.patient,condition=condition,library=library,cutoff=cutoff,
                    label=group[0],category=group[1],n_cells=len(take),
                    median_nCount=float(take.nCount.median()),median_nFeature=float(take.nFeature.median()),
                    median_pct_mt_visible=float(take.pct_mt_visible.median())))
            errors.append(dict(sample=s.sample,patient=s.patient,condition=condition,library=library,cutoff=cutoff,
                n_cells=len(p),n_marker_known=int((~undecided).sum()),n_marker_wrong=int(((~undecided)&(seed!=truth.L1)).sum()),
                n_new_correct=int((category=='new_correct').sum()),n_new_incorrect=int((category=='new_incorrect').sum()),
                n_terminal_unknown=int(unknown.sum())))
            # Within-original-class QC differences prevent class composition alone
            # from being mistaken for an Unknown quality signal.
            for label,take in data.groupby('truth'):
                if not take.unknown.any() or take.unknown.all():continue
                a=take[take.unknown];b=take[~take.unknown]
                rows.append(dict(sample=s.sample,patient=s.patient,condition=condition,library=library,cutoff=cutoff,
                    label=label,category='Unknown_minus_called_within_class',n_cells=len(take),
                    median_nCount=float(a.nCount.median()-b.nCount.median()),
                    median_nFeature=float(a.nFeature.median()-b.nFeature.median()),
                    median_pct_mt_visible=float(a.pct_mt_visible.median()-b.pct_mt_visible.median())))
    table=pd.DataFrame(rows);error=pd.DataFrame(errors)
    table.to_csv(dest/'QC_by_sample_truth_and_call_status.csv',index=False)
    error.to_csv(dest/'retained_seed_and_DL_new_errors.csv',index=False)
    contrasts=table[table.category=='Unknown_minus_called_within_class']
    patient=contrasts.groupby(['patient','condition'])[['median_nCount','median_nFeature','median_pct_mt_visible']].mean().reset_index()
    patient.to_csv(dest/'within_truth_class_QC_patient_differences.csv',index=False)
    fig,axs=plt.subplots(1,3,figsize=(13,4),layout='constrained')
    for ax,key,label in zip(axs,['median_nCount','median_nFeature','median_pct_mt_visible'],['Counts','Detected genes','Visible mitochondrial %']):
        for i,condition in enumerate(['fixed_glioma','training_patient_selected']):
            v=patient.loc[patient.condition==condition,key].to_numpy()
            if len(v):ax.scatter(np.full(len(v),i)+np.linspace(-.09,.09,len(v)),v,s=12,alpha=.7,color=['#0072B2','#D55E00'][i])
        ax.axhline(0,color='#555',lw=1);ax.set_xticks([0,1],['Fixed marker','Selected marker'],rotation=15)
        ax.set(title=label,ylabel='Unknown minus called; class-matched median difference')
    fig.suptitle('QC association, not a diagnosis: each point is a patient\nOnly samples/classes containing both called and Unknown cells contribute')
    for ext in ['png','pdf']:fig.savefig(dest/f'Unknown_QC_association.{ext}',dpi=200,bbox_inches='tight')
    plt.close(fig)
    (dest/'INTERPRETATION.md').write_text('Unknown is an abstention status. The QC comparisons condition on original author L1 class and average within patients. They are descriptive associations, not proof of doublets, poor-quality cells or a novel cell type. No doublet scores were available in this audited QC export. Mitochondrial percentage refers only to mitochondrial genes visible in the supplied matrix. Wrong known marker seeds are retained by the original algorithm; newly filled wrong predictions are counted separately.\n')
    write_json(dest/'manifest.json',dict(status='completed',n_primary_samples=97,sources=sources,
        no_new_fit=True,no_doublet_or_novel_type_claim=True,
        files={p.name:sha(p) for p in dest.iterdir() if p.suffix in ['.csv','.png','.pdf','.md']},
        job=os.environ['SLURM_JOB_ID'],source_sha256=sha(__file__),completed_at=utc()))
    complete(dest)

if __name__=='__main__':run()
