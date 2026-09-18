"""Patient-level class performance, fixed-display confusion and actual DL changes."""
import json
import os
from common import OUT, L1, ROUTES, require_slurm, checked, sha, write_json, utc

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d=OUT/'summary';assert checked(d,'aggregate_manifest.json','AGGREGATE_COMPLETE')
    cohort=pd.read_csv(OUT/'protocol/cohort.csv')
    frames=[]
    for sample in cohort['sample']:
        p=OUT/'GBM'/sample/'hvg2000/evaluation/per_class.csv.gz'
        f=pd.read_csv(p,dtype={'cutoff':str})
        frames.append(f[(f.library=='CM2_glioma_other')&(f.cutoff=='mean')&(f.stage=='terminal090')&(f.family=='native_R_budget')])
    per=pd.concat(frames,ignore_index=True)
    # A class absent from a patient's samples is not a patient-level observation.
    patient=per[per.support>0].groupby(['patient','primary','route','label'])[['F1','precision','recall','support']].mean().reset_index()
    patient.to_csv(d/'original2000_fixed_marker_per_class_patient.csv',index=False)
    mean=patient[patient.primary].groupby(['route','label'])[['F1','precision','recall']].mean().reset_index()
    mean.to_csv(d/'original2000_fixed_marker_per_class_mean.csv',index=False)
    f1=mean.pivot(index='label',columns='route',values='F1').reindex(index=L1,columns=ROUTES)
    fig,ax=plt.subplots(figsize=(8.5,5.8),layout='constrained')
    im=ax.imshow(f1.to_numpy(),cmap='viridis',vmin=0,vmax=1,aspect='auto')
    ax.set_xticks(range(4),['PCA / SNN','PCA / HDBSCAN','UMAP / SNN','UMAP / HDBSCAN'],rotation=20,ha='right')
    ax.set_yticks(range(11),L1);ax.set_title('HVG2000: per-class terminal F1\nFixed glioma marker; primary-cohort patient means')
    for i in range(11):
        for j in range(4):
            v=f1.iloc[i,j]
            ax.text(j,i,'NA' if pd.isna(v) else f'{v:.2f}',ha='center',va='center',fontsize=9,color='white' if pd.isna(v) or v<.55 else 'black')
    fig.colorbar(im,ax=ax,label='F1')
    for ext in ['png','pdf']:fig.savefig(d/f'per_class_original2000.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)
    display=json.loads((OUT/'protocol/input_audit.json').read_text())['display_sample']
    cm=pd.read_csv(OUT/'GBM'/display/'hvg2000/evaluation/confusions.csv.gz',dtype={'cutoff':str})
    cm=cm[(cm.library=='CM2_glioma_other')&(cm.cutoff=='mean')&(cm.stage=='terminal090')&(cm.family=='native_R_budget')]
    columns=L1+sorted(set(cm.prediction)-set(L1))
    fig,axs=plt.subplots(2,2,figsize=(16,12),layout='constrained')
    for ax,route in zip(axs.flat,ROUTES):
        t=cm[cm.route==route].pivot_table(index='truth',columns='prediction',values='n',aggfunc='sum',fill_value=0).reindex(index=L1,columns=columns,fill_value=0)
        assert t.to_numpy().sum()==cm.loc[cm.route==route,'n'].sum()
        den=t.sum(axis=1).replace(0,np.nan);v=t.div(den,axis=0)
        im=ax.imshow(v.to_numpy(),cmap='Blues',vmin=0,vmax=1,aspect='auto')
        ax.set_xticks(range(len(columns)),columns,rotation=70,ha='right',fontsize=7)
        ax.set_yticks(range(len(L1)),L1,fontsize=8);ax.set_title(route)
        ax.set(xlabel='Final prediction; Unknown and unmapped retained',ylabel='Original author label')
    fig.suptitle(f'{display}: fixed HVG2000 glioma-marker terminal confusion\nRows normalized by all truth-class cells; missing truth classes shown blank')
    fig.colorbar(im,ax=list(axs.flat),label='Fraction of truth-class cells',shrink=.6)
    for ext in ['png','pdf']:fig.savefig(d/f'display_sample_confusions.{ext}',dpi=200,bbox_inches='tight')
    plt.close(fig)
    a=pd.read_csv(d/'all_annotation_metrics.csv.gz',dtype={'cutoff':str})
    a=a[(a.stage=='terminal090')&(a.library=='CM2_glioma_other')&(a.cutoff=='mean')]
    columns=['sample','patient','primary','budget','route','dl_status','training_executed','n_cells','n_known','n_pool',
        'n_training_classes','n_new_correct','n_new_incorrect','n_initially_wrong_retained','coverage','unknown_rate']
    a[columns].to_csv(d/'fixed_marker_DL_change_audit.csv',index=False)
    write_json(d/'diagnostics_manifest.json',dict(status='completed',display_sample=display,
        missing_classes_remain_errors=True,known_marker_labels_retained_by_design=True,
        propagation='Original DL only fills Undecided; wrong known seeds are retained. New wrong fills are reported explicitly.',
        files={p.name:sha(p) for p in d.iterdir() if p.stem in ['per_class_original2000','display_sample_confusions','fixed_marker_DL_change_audit','original2000_fixed_marker_per_class_patient','original2000_fixed_marker_per_class_mean']},
        job=os.environ['SLURM_JOB_ID'],source_sha256=sha(__file__),completed_at=utc()))

if __name__=='__main__':run()
