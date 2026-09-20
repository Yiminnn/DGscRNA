"""Author truth is read only after the complete candidate prediction set is frozen."""
from pathlib import Path
import json
import os
import sys
from a2_common import CODE, OUT, REFERENCE, L1, verify_source, verify_terminal, sha, utc, write_json, checked, complete, require_slurm


def run(directory):
    require_slurm();protocol=verify_source()
    import numpy as np
    import pandas as pd
    from sklearn.metrics import precision_recall_fscore_support
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    directory=Path(directory);cfg=json.loads((directory/'config.json').read_text())
    assert checked(directory,'fit_manifest.json','FIT_COMPLETE')
    fm=json.loads((directory/'fit_manifest.json').read_text())
    assert fm['input_signature']==cfg['input_signature']
    route=directory/'cellwise_seed';dest=directory/'evaluation';dest.mkdir(exist_ok=True)
    signatures={arm['id']:sha(route/'terminal'/arm['id']/'terminal_manifest.json') for arm in cfg['arms']}
    assert signatures==fm['terminal_manifests']
    if checked(dest):
        old=json.loads((dest/'manifest.json').read_text())
        assert old['terminal_manifest_hashes']==signatures and old['source_bundle_sha256']==sha(CODE/'SOURCE_MANIFEST.json')
        for name,digest in old['outputs'].items():assert sha(dest/name)==digest,name
        return
    sample,budget=cfg['sample'],cfg['budget']
    im=json.loads((REFERENCE/'inputs'/sample/'input_manifest.json').read_text())
    truthpath=REFERENCE/'evaluation_inputs'/sample/'truth.csv.gz'
    assert sha(truthpath)==im['evaluation_files']['truth.csv.gz']==cfg['input']['truth_sha256']
    truth=pd.read_csv(truthpath,dtype=str,keep_default_na=False)
    cells=pd.read_csv(route/'cells.csv',dtype=str,keep_default_na=False)
    assert np.array_equal(cells.cell_id,truth.cell_id)
    mpath=REFERENCE/'markers/panel_L1_mapping.csv'
    assert sha(mpath)==protocol['mapping_sha256']
    mapping=pd.read_csv(mpath,dtype=str,keep_default_na=False)
    mapping=mapping[mapping.library=='CM2_glioma_other'];lookup=dict(zip(mapping.panel,mapping.L1))
    y=truth.L1.to_numpy(dtype=str)
    present_labels=sorted(set(y));rows=[];perclass=[];confusions=[];pictures=[]
    coarse=lambda a:np.asarray(['Neuron' if v in ['Excitatory neuron','Inhibitory neuron','AMBIGUOUS_NEURON'] else v for v in a])
    coarse_labels=[v for v in L1 if v not in ['Excitatory neuron','Inhibitory neuron']]+['Neuron'];yc=coarse(y)
    for arm in cfg['arms']:
        aid=arm['id'];tm=verify_terminal(route,aid)
        z=np.load(route/'terminal'/aid/'terminal.npz',allow_pickle=False)
        assert len(z['initial'])==len(y)
        known=z['initial']!='Undecided'
        assert np.array_equal(z['initial'][known],z['final090'][known])
        stage_preds={}
        for stage,key in [('marker_only','initial'),('terminal090','final090'),('terminal070','final070')]:
            native=z[key]
            pred=np.asarray([lookup.get(v,'Unknown' if v in ['Unknown','Undecided','Noise',''] else 'UNMAPPABLE') for v in native])
            precision,recall,f1,support=precision_recall_fscore_support(y,pred,labels=L1,zero_division=0)
            cp,cr,cf,cs=precision_recall_fscore_support(yc,coarse(pred),labels=coarse_labels,zero_division=0)
            abstain=np.isin(native,['Unknown','Undecided','Noise','']);mapped=np.isin(pred,L1)
            context=dict(sample=sample,patient=im['patient'],primary=im['primary'],budget=budget,
                route='cellwise_seed',library='CM2_glioma_other',arm_id=aid,lambda_value=arm['lambda'],stage=stage,
                family='cluster_DEG_seed_replacement',n_cells=len(y),seed=42,DL_features=tm['DL_features'])
            result=dict(**context,status='completed',dl_status=tm['dl_status'],training_executed=tm['training_executed'],
                macroF1_present=float(f1[support>0].mean()),macroF1_fixed11=float(f1.mean()),
                weightedF1=float((f1*support).sum()/support.sum()),accuracy=float((pred==y).mean()),
                coverage=float((~abstain).mean()),unknown_rate=float(abstain.mean()),mapped_coverage=float(mapped.mean()),
                off_vocabulary_rate=float((~mapped & ~abstain).mean()),
                coarse10_macroF1_present=float(cf[cs>0].mean()),coarse10_macroF1_fixed10=float(cf.mean()),
                n_known=tm['n_known'],n_pool=tm['n_pool'],n_training_classes=tm['n_training_classes'],
                n_new_correct=int((~known & (pred==y)).sum()),n_new_incorrect=int((~known & ~abstain & (pred!=y)).sum()),
                n_initially_wrong_retained=int((known & (pred!=y)).sum()))
            rows.append(result);stage_preds[stage]=(pred,result)
            for index,label in enumerate(L1):
                perclass.append(dict(**context,label=label,precision=float(precision[index]),recall=float(recall[index]),F1=float(f1[index]),support=int(support[index])))
            for record in pd.DataFrame({'truth':y,'prediction':pred}).value_counts().reset_index(name='n').to_dict('records'):
                confusions.append(dict(**context,**record))
        pictures.append((arm,tm,stage_preds))
        z.close()
    pd.DataFrame(rows).to_csv(dest/'metrics.csv',index=False)
    pd.DataFrame(perclass).to_csv(dest/'per_class.csv.gz',index=False,compression='gzip')
    pd.DataFrame(confusions).to_csv(dest/'confusions.csv.gz',index=False,compression='gzip')

    display=REFERENCE/'GBM'/sample/'hvg2000/UMAP2.csv'
    assert sha(display)==cfg['input']['display_sha256']
    xy=pd.read_csv(display,index_col=0);assert list(xy.index)==list(truth.cell_id)
    xy=xy.to_numpy();cm=plt.get_cmap('tab20');colors={label:cm(i) for i,label in enumerate(L1)}
    colors.update(Unknown='#bdbdbd',UNMAPPABLE='#7b614c',AMBIGUOUS_NEURON='#8c8c33',NO_L1_COUNTERPART='#7b614c')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42})
    fig,axes=plt.subplots(2,6,figsize=(21,7.5),layout='constrained')
    shown=set(y)
    def panel(ax,pred,title):
        pred=np.asarray(pred);shown.update(pred)
        for label in sorted(set(pred)):
            take=pred==label;ax.scatter(xy[take,0],xy[take,1],s=2.5,c=[colors.get(label,'#7b614c')],linewidths=0,rasterized=True,alpha=.8)
        ax.set_title(title,fontsize=8);ax.set_xticks([]);ax.set_yticks([])
        ax.set_xlabel('Shared original HVG2000 UMAP1',fontsize=6);ax.set_ylabel('UMAP2',fontsize=6)
        for spine in ax.spines.values():spine.set_visible(False)
    panel(axes[0,0],y,'Original author L1');panel(axes[1,0],y,'Original author L1')
    for column,(arm,tm,stages) in enumerate(pictures,1):
        for row,stage in enumerate(['marker_only','terminal090']):
            pred,metric=stages[stage]
            descriptor='Cell-wise seed' if row==0 else 'Terminal 0.90: '+tm['dl_status']
            panel(axes[row,column],pred,f"lambda {arm['lambda']:g} | {descriptor}\nmacro-F1 {metric['macroF1_present']:.3f}; coverage {metric['coverage']:.1%}")
    fig.suptitle(f'{sample} | {budget} | CM2_glioma_other\nCell-wise expression seeds replace cluster-DEG seeds; displayed UMAP is not a model input',fontsize=12)
    fig.legend(handles=[Line2D([],[],marker='o',color='none',markerfacecolor=colors.get(label,'#7b614c'),markersize=5,label=label) for label in sorted(shown)],
               loc='outside lower center',ncol=7,fontsize=8,frameon=False)
    for suffix in ['png','pdf']:fig.savefig(dest/f'candidate_seed_terminal_grid.{suffix}',dpi=180,bbox_inches='tight')
    plt.close(fig)
    outputs=['metrics.csv','per_class.csv.gz','confusions.csv.gz','candidate_seed_terminal_grid.png','candidate_seed_terminal_grid.pdf']
    write_json(dest/'manifest.json',dict(status='completed',sample=sample,budget=budget,n_cells=len(y),n_candidates=5,
        n_metric_rows=len(rows),terminal_manifest_hashes=signatures,truth_sha256=sha(truthpath),mapping_sha256=sha(mpath),
        display_sha256=sha(display),every_candidate_plotted=True,all_cell_denominator=True,
        source_bundle_sha256=sha(CODE/'SOURCE_MANIFEST.json'),no_clustering_metric_fabricated=True,
        outputs={name:sha(dest/name) for name in outputs},job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest);print('A2_EVALUATION_COMPLETE',sample,budget,len(rows),flush=True)


if __name__=='__main__':run(sys.argv[1])
