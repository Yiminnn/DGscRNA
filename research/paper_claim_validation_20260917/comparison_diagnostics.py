"""Selected comparator class errors and fixed-coordinate maps, all cells retained."""
import json
import os
from common import OUT, L1, require_slurm, checked, complete, sha, write_json, utc

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from prediction_helpers import read_native,map_labels
    from sklearn.metrics import precision_recall_fscore_support

    root=OUT/'comparison_summary';assert checked(root)
    dest=root/'diagnostics';dest.mkdir(exist_ok=True)
    cohort=pd.read_csv(OUT/'protocol/cohort.csv')
    cv=pd.read_csv(root/'patient_heldout_results.csv',dtype={'cutoff':str})
    methods=['DG-scRNA','scType','scCATCH','SCINA','SingleR','scDeepSort']
    allclasses=[];allconfusions=[]
    for sample in cohort.itertuples():
        for method in methods:
            path=OUT/'GBM'/sample.sample/'hvg2000/evaluation' if method=='DG-scRNA' else OUT/'comparators'/method/sample.sample/'evaluation'
            classes=pd.read_csv(path/'per_class.csv.gz',dtype={'cutoff':str})
            confusion=pd.read_csv(path/'confusions.csv.gz',dtype={'cutoff':str})
            chosen=cv[(cv.patient==sample.patient)&(cv.method==method)]
            for config in chosen.to_dict('records'):
                if config['cohort']=='primary97' and not sample.primary:continue
                def select(frame):
                    mask=np.ones(len(frame),dtype=bool)
                    for key in ['budget','route','library','cutoff']:mask &= frame[key].eq(config[key]).to_numpy()
                    if method=='DG-scRNA':mask &= frame.stage.eq('terminal090').to_numpy()
                    take=frame.loc[mask].copy();assert len(take),(sample.sample,method,config)
                    take['cohort']=config['cohort'];take['method']=method
                    return take
                c=select(classes);f=select(confusion)
                assert len(c)==len(L1) and c.support.sum()==sample.n_cells and f.n.sum()==sample.n_cells
                allclasses.append(c);allconfusions.append(f)
    classes=pd.concat(allclasses,ignore_index=True);confusion=pd.concat(allconfusions,ignore_index=True)
    classes.to_csv(dest/'selected_per_sample_classes.csv.gz',index=False,compression='gzip')
    confusion.to_csv(dest/'selected_all_cell_confusions.csv.gz',index=False,compression='gzip')
    patient=classes[classes.support>0].groupby(['cohort','method','patient','label'])[['precision','recall','F1']].mean().reset_index()
    patient.to_csv(dest/'selected_per_patient_classes.csv',index=False)
    mean=patient.groupby(['cohort','method','label'])[['precision','recall','F1']].mean().reset_index()
    mean.to_csv(dest/'selected_class_means.csv',index=False)
    grid=mean[mean.cohort=='primary97'].pivot(index='label',columns='method',values='F1').reindex(index=L1,columns=methods)
    fig,ax=plt.subplots(figsize=(10,6),layout='constrained');im=ax.imshow(grid,cmap='viridis',vmin=0,vmax=1,aspect='auto')
    ax.set_xticks(range(len(methods)),methods,rotation=25,ha='right');ax.set_yticks(range(len(L1)),L1)
    for i in range(len(L1)):
        for j in range(len(methods)):
            v=grid.iloc[i,j]
            ax.text(j,i,'NA' if pd.isna(v) else f'{v:.2f}',ha='center',va='center',fontsize=8,color='white' if pd.isna(v) or v<.55 else 'black')
    ax.set_title('Heldout-patient per-class F1: every cell retained\nMarker methods | labelled-reference SingleR | atlas GNN scDeepSort')
    fig.colorbar(im,ax=ax,label='F1')
    for ext in ['png','pdf']:fig.savefig(dest/f'heldout_class_F1.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)
    sample=json.loads((OUT/'protocol/input_audit.json').read_text())['display_sample']
    sid=cohort[cohort['sample']==sample].iloc[0]
    truth=pd.read_csv(OUT/'evaluation_inputs'/sample/'truth.csv.gz',dtype=str,keep_default_na=False)
    coords=pd.read_csv(OUT/'GBM'/sample/'hvg2000/UMAP2.csv',index_col=0)
    assert list(coords.index)==list(truth.cell_id)
    chosen=cv[(cv.cohort=='primary97')&(cv.patient==sid.patient)]
    panels=[('Author L1',truth.L1.to_numpy())];proof=[]
    for method in methods:
        cfg=chosen[chosen.method==method].iloc[0]
        pred,source=read_native(method,sample,cfg.budget,cfg.route,cfg.library,cfg.cutoff)
        assert np.array_equal(pred.cell_id,truth.cell_id)
        p=map_labels(method,cfg.library,pred.prediction)
        _,_,f1,su=precision_recall_fscore_support(truth.L1,p,labels=L1,zero_division=0)
        # Numerical parity against existing all-cell evaluation, before visualization.
        src=OUT/'GBM'/sample/'hvg2000/evaluation' if method=='DG-scRNA' else OUT/'comparators'/method/sample/'evaluation'
        saved=pd.read_csv(src/'metrics.csv',dtype={'cutoff':str})
        keep=np.ones(len(saved),dtype=bool)
        for key in ['budget','route','library','cutoff']:keep &= saved[key].eq(cfg[key]).to_numpy()
        if method=='DG-scRNA':keep &= saved.stage.eq('terminal090').to_numpy()
        row=saved[keep];assert len(row)==1
        assert np.isclose(float(f1[su>0].mean()),row.iloc[0].macroF1_present,rtol=0,atol=1e-14)
        proof.append(dict(method=method,source=str(source),sha256=sha(source),mapping_parity=True))
        panels.append((method+'\n'+cfg.library+' / '+str(cfg.cutoff),p))
    colors={v:plt.get_cmap('tab20')(i) for i,v in enumerate(L1)}
    colors.update(Unknown='#bdbdbd',UNMAPPABLE='#7b614c',AMBIGUOUS_NEURON='#8c8c33',NO_L1_COUNTERPART='#7b614c')
    fig,axs=plt.subplots(2,4,figsize=(18,9),layout='constrained')
    for ax,(title,labs) in zip(axs.flat,panels):
        for lab in sorted(set(labs)):
            mask=labs==lab
            ax.scatter(coords.iloc[mask,0],coords.iloc[mask,1],s=2,c=[colors.get(lab,'#7b614c')],linewidths=0,rasterized=True)
        ax.set(title=title,xlabel='Fixed display UMAP1',ylabel='Fixed display UMAP2');ax.set_xticks([]);ax.set_yticks([])
    axs.flat[-1].axis('off')
    fig.legend(handles=[Line2D([],[],marker='o',color='none',markerfacecolor=colors[v],label=v,markersize=5) for v in L1+['Unknown','UNMAPPABLE','AMBIGUOUS_NEURON']],
               loc='outside lower center',ncol=7,fontsize=8,frameon=False)
    fig.suptitle(f'{sample}: marker/threshold choices made on training patients\nIdentical display coordinates; SingleR and scDeepSort have distinct reference information')
    for ext in ['png','pdf']:fig.savefig(dest/f'display_sample_methods.{ext}',dpi=200,bbox_inches='tight')
    plt.close(fig)
    summary=pd.read_csv(root/'patient_heldout_summary.csv')
    primary=summary[summary.cohort=='primary97'].set_index('method').loc[methods]
    stats=pd.read_csv(root/'paired_patient_comparisons.csv');stats=stats[stats.cohort=='primary97'].set_index('method').loc[methods[1:]]
    fig,axs=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    axs[0].bar(range(6),primary.macroF1_present_mean,color=['#0072B2']*4+['#D55E00','#009E73'])
    axs[0].plot(range(6),primary.coverage_mean,'ko--',ms=4,label='Native-call coverage')
    axs[0].set_xticks(range(6),methods,rotation=30,ha='right');axs[0].set(ylim=(0,1),ylabel='Patient mean',title='Terminal F1 and call coverage');axs[0].legend(frameon=False)
    axs[1].errorbar(stats.mean_delta,range(5),xerr=np.vstack([stats.mean_delta-stats.CI95_low,stats.CI95_high-stats.mean_delta]),fmt='o',color='#0072B2')
    axs[1].axvline(0,color='#555',ls='--');axs[1].set_yticks(range(5),methods[1:]);axs[1].set(xlabel='Competitor minus DG-scRNA macro-F1',title='Paired patients: conditional 95% bootstrap interval')
    fig.suptitle('Primary97 / 55 patients: retrospective label-heldout comparisons\nMatched marker opportunities; reference-information conditions shown explicitly')
    for ext in ['png','pdf']:fig.savefig(dest/f'heldout_method_comparison.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)
    write_json(dest/'manifest.json',dict(status='completed',display_sample=sample,mapping_parity=proof,
        all_cells_in_confusions=True,files={p.name:sha(p) for p in dest.iterdir() if p.suffix in ['.csv','.gz','.png','.pdf']},
        job=os.environ['SLURM_JOB_ID'],source_sha256=sha(__file__),completed_at=utc()))
    complete(dest)

if __name__=='__main__':run()
