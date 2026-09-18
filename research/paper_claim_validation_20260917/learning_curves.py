"""Plot saved native-MLP histories and epoch endpoints, without new fitting."""
import json
import os
from pathlib import Path
from common import OUT,ROUTES,require_slurm,checked,sha,write_json,complete,utc


def run(dest=None):
    require_slurm()
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from dl_controls import CONFIGS,LIBRARIES
    dest=Path(dest or OUT/'controls_summary');dest.mkdir(exist_ok=True)
    if checked(dest,'learning_curves_manifest.json','LEARNING_CURVES_COMPLETE'):return
    pilots=(OUT/'protocol/pilot_samples.txt').read_text().split()
    histories=[];states=[];endpoints=[];sources={}
    keys=['sample','budget','route','library','control']
    for sample in pilots:
        for budget in ['hvg2000','hvg5000']:
            for control in CONFIGS:
                folder=OUT/'GBM_DL_controls'/sample/budget/control
                assert checked(folder)
                metric=pd.read_csv(folder/'metrics.csv',dtype={'cutoff':str})
                for route in ROUTES:
                    sm=json.loads((folder/route/'score_manifest.json').read_text())
                    arms={aid:arm for aid,arm in sm['arms'].items()
                          if arm['library'] in LIBRARIES and arm['cutoff']=='mean'}
                    assert len(arms)==2
                    for aid,arm in arms.items():
                        p=folder/route/'terminal'/aid
                        assert checked(p,'terminal_manifest.json','TERMINAL_COMPLETE')
                        tm=json.loads((p/'terminal_manifest.json').read_text())
                        tr=json.loads((p/'training_manifest.json').read_text())
                        assert sha(p/'training_manifest.json')==tm['training_manifest_sha256']
                        key=dict(sample=sample,budget=budget,route=route,library=arm['library'],control=control)
                        params=tr['params'];assert all(params[k]==v for k,v in CONFIGS[control].items())
                        state=dict(**key,epochs=params['epochs'],model_seed=params['model_seed'],
                            architecture=str(params['architecture']),training_executed=tr['training_executed'],
                            dl_status=tr['dl_status'],n_known=tr['n_known'],n_pool=tr['n_pool'],
                            n_training_classes=tr['n_training_classes'],n_train=tr.get('n_train',0),
                            n_validation=tr.get('n_validation',0),
                            known_seed_validation_accuracy=tr.get('known_seed_validation_accuracy'))
                        states.append(state)
                        sources[str((p/'training_manifest.json').relative_to(OUT))]=sha(p/'training_manifest.json')
                        if tr['training_executed']:
                            history_path=p/'training_history.json'
                            assert sha(history_path)==tr['outputs']['training_history.json']
                            history=json.loads(history_path.read_text())
                            assert [h['epoch'] for h in history]==list(range(1,params['epochs']+1))
                            assert all(h['n']==tr['n_train'] and np.isfinite(h['mean_loss']) and
                                       0<=h['accuracy']<=1 for h in history)
                            histories.extend(dict(**key,**h) for h in history)
                        if control in ['epochs5','original','epochs20']:
                            m=metric[(metric.route==route)&(metric.library==arm['library'])&metric.stage.eq('final090')]
                            assert len(m)==1
                            endpoints.append(dict(**state,macroF1_present=float(m.iloc[0].macroF1_present),
                                                  coverage=float(m.iloc[0].coverage)))
    history=pd.DataFrame(histories);status=pd.DataFrame(states);end=pd.DataFrame(endpoints)
    assert len(status)==3*2*9*4*2 and len(end)==3*2*3*4*2
    history.to_csv(dest/'MLP_learning_history.csv',index=False)
    status.to_csv(dest/'MLP_learning_training_status.csv',index=False)
    end.to_csv(dest/'MLP_epoch_endpoint_checks.csv',index=False)
    checks=[]
    for key,g in status.groupby(keys[:-1],sort=True):
        a=history
        for field,value in zip(keys[:-1],key):a=a[a[field]==value]
        epoch_runs=g[g.control.isin(['epochs5','original','epochs20'])]
        assert epoch_runs.training_executed.nunique()==1
        if epoch_runs.training_executed.iloc[0]:
            longest=a[a.control=='epochs20'].set_index('epoch')
            for label,n in [('epochs5',5),('original',10)]:
                shorter=a[a.control==label].set_index('epoch')
                for field in ['n','mean_loss','accuracy','legacy_batch_size_weighted_sum']:
                    np.testing.assert_allclose(shorter[field],longest.loc[range(1,n+1),field],rtol=1e-10,atol=1e-12)
            checks.append(dict(zip(keys[:-1],key),status='saved_5_10_20_epoch_prefixes_match'))
        else:
            checks.append(dict(zip(keys[:-1],key),status='no_training_curve',dl_status=epoch_runs.iloc[0].dl_status))
    pd.DataFrame(checks).to_csv(dest/'MLP_epoch_history_prefix_parity.csv',index=False)
    colors=['#0072B2','#D55E00','#009E73','#CC79A7']
    styles=['-','--'];handles=[]
    for route,color in zip(ROUTES,colors):
        for lib,style in zip(LIBRARIES,styles):
            handles.append(Line2D([],[],color=color,linestyle=style,label=route+' / '+lib,linewidth=1.2))
    for kind,title,ylabel,name in [
        ('loss','Saved training loss: original width and initialization, up to 20 epochs',
         'Mean legacy Softmax + cross-entropy loss','MLP_training_learning_curves'),
        ('endpoint','Annotation after 5, 10 and 20 epochs: original width and initialization',
         'All-cell author-L1 terminal macro-F1','MLP_epoch_endpoints')]:
        fig,axes=plt.subplots(3,2,figsize=(13,11),layout='constrained')
        for i,sample in enumerate(pilots):
            for j,budget in enumerate(['hvg2000','hvg5000']):
                ax=axes[i,j];plotted=0
                for route,color in zip(ROUTES,colors):
                    for lib,style in zip(LIBRARIES,styles):
                        source=history if kind=='loss' else end
                        g=source[(source['sample']==sample)&(source.budget==budget)&
                                 (source.route==route)&(source.library==lib)]
                        if kind=='loss':
                            g=g[g.control=='epochs20'].sort_values('epoch')
                            if len(g):ax.plot(g.epoch,g.mean_loss,color=color,linestyle=style,linewidth=1.2);plotted+=1
                        else:
                            g=g.sort_values('epochs')
                            assert len(g)==3
                            ax.plot(g.epochs,g.macroF1_present,'o',color=color,linestyle=style,ms=3,linewidth=1.2)
                if kind=='loss' and not plotted:ax.text(.5,.5,'No trainable refinement in these conditions',ha='center',va='center',transform=ax.transAxes)
                ax.axvline(10,color='#777',linestyle=':',linewidth=.8)
                ax.set(title=sample+' / '+budget,xlabel='Epoch',ylabel=ylabel)
                if kind=='endpoint':ax.set_xticks([5,10,20]);ax.set_ylim(0,1)
                ax.spines[['top','right']].set_visible(False)
        fig.suptitle(title+'\nThree count-selected pilots; fixed clusters, marker seeds and split42')
        fig.legend(handles=handles,loc='outside lower center',ncol=2,fontsize=7,frameon=False)
        for ext in ['png','pdf']:fig.savefig(dest/(name+'.'+ext),dpi=200,bbox_inches='tight')
        plt.close(fig)
    note='''# Saved learning curves

No model is refitted for this diagnostic. Every history and terminal state from
the 54 completed MLP-control units is exported, including unavailable/no-op
conditions. Epoch-specific losses are the actual per-observation weighted training
loss, not the historical sum that multiplied each batch by256. The historical
Softmax-to-CrossEntropyLoss combination remains unchanged.

The main loss figure uses the20-epoch arm at original width and model seed42.
The saved5-epoch and10-epoch histories are checked against its first5/10epochs.
No-op or structurally untrainable conditions have no invented loss curve. The
companion plot shows all-cell terminal annotation F1 at5/10/20epochs for every
route and both frozen marker contexts, including valid no-op endpoints.

Training targets are marker pseudo-labels. Final known-seed validation accuracy is
exported per run; per-epoch validation losses were not recorded and are not
fabricated. Agreement with seed labels does not validate biological correctness,
and a decreasing training loss does not establish optimal stopping or convergence
on true cell types. These are three prespecified size pilots, not55 independent
patients or a new all-cohort epoch-selection experiment.
'''
    (dest/'MLP_LEARNING_CURVE_INTERPRETATION.md').write_text(note)
    names=['MLP_learning_history.csv','MLP_learning_training_status.csv','MLP_epoch_endpoint_checks.csv',
           'MLP_epoch_history_prefix_parity.csv','MLP_LEARNING_CURVE_INTERPRETATION.md']
    names += [base+'.'+ext for base in ['MLP_training_learning_curves','MLP_epoch_endpoints'] for ext in ['png','pdf']]
    write_json(dest/'learning_curves_manifest.json',dict(status='completed',n_conditions=len(status),
        n_histories=int(status.training_executed.sum()),n_epoch_rows=len(history),n_epoch_endpoint_conditions=len(end),
        no_new_fits=True,scope='Three count-selected size pilots; marker-pseudo-label training',
        source_manifests=sources,files={name:sha(dest/name) for name in names},
        source_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest,'learning_curves_manifest.json','LEARNING_CURVES_COMPLETE')


if __name__=='__main__':run()
