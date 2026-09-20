"""Verify the full frozen B_GBM scope and summarize genuine validation/neighbor results."""
from pathlib import Path
import json,os,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from run_controls import OUT,OLD,OLD_CODE,LIBRARIES,PILOTS,BUDGETS,tasks,sha,utc,write_json,checked,complete,require_slurm

def refresh_status():
    done=0;failed=[]
    for i,cfg in enumerate(tasks()):
        if cfg['task']=='neighbors':
            p=OUT/'neighbors'/cfg['sample']/cfg['budget']/cfg['name'];ok=checked(p)
        else:
            p=OUT/'learning'/cfg['sample']/cfg['budget']/f"seed{cfg['model_seed']}"/LIBRARIES[0]
            ok=checked(p,'manifest.json','LEARNING_COMPLETE')
            if ok and not json.loads((p/'manifest.json').read_text())['training_executed']:
                ok=checked(p.parent/LIBRARIES[1],'manifest.json','LEARNING_COMPLETE')
        done+=int(ok)
        f=OUT/'failures'/f'task_{i}.json'
        if f.exists() and not ok:failed.append(str(f))
    p=OUT/'status.json';previous=json.loads(p.read_text()) if p.exists() else {}
    previous.update(stage='B_GBM',status='compute_complete_pending_aggregation' if done==84 else ('running_with_failures' if failed else 'running'),updated_at=utc(),completed=done,remaining=84-done,
        failures=failed,evidence=[str(OUT/'tasks.json'),str(OUT/'verification/parity.json')])
    write_json(p,previous);return previous

def aggregate():
    require_slurm()
    import numpy as np,pandas as pd,matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    state=refresh_status();assert state['completed']==84,state
    assert checked(OUT/'verification','parity.json','PARITY_COMPLETE')
    parity_proof=json.loads((OUT/'verification/parity.json').read_text())
    # Post-pilot orchestration/figure improvements do not replace the numerical source gate.
    for name in ['prepare_neighbors.R','score_neighbors.R','refine_learning.py']:
        assert sha(HERE/name)==parity_proof['source_hashes'][name]
    deriv=json.loads((HERE/'SOURCE_DERIVATION.json').read_text())
    for name,digest in deriv.items():assert sha(OLD_CODE/name)==digest
    neighbors=[];metrics=[];histories=[];learning_states=[];proof=[];artifacts={}
    for cfg in tasks():
        if cfg['task']=='neighbors':
            p=OUT/'neighbors'/cfg['sample']/cfg['budget']/cfg['name'];m=json.loads((p/'manifest.json').read_text())
            for name,digest in m['files'].items():assert sha(p/name)==digest
            score=p/(cfg['space']+'_'+cfg['method'])
            assert checked(score,'score_manifest.json','SCORE_COMPLETE')
            sm=json.loads((score/'score_manifest.json').read_text())
            assert sha(score/'initial_calls.csv.gz')==sm['initial_sha256']
            for aid,a in sm['arms'].items():
                if a['library'] not in LIBRARIES or a['cutoff']!='mean':continue
                td=score/'terminal'/aid
                assert checked(td,'terminal_manifest.json','TERMINAL_COMPLETE')
                tm=json.loads((td/'terminal_manifest.json').read_text())
                assert sha(td/'training_manifest.json')==tm['training_manifest_sha256']
                assert sha(td/'terminal.npz')==tm['terminal_sha256']
                assert sha(td/'predictions.csv.gz')==tm['predictions_sha256']
            z=pd.read_csv(p/'metrics.csv');assert len(z)==6
            assert set(z.stage)=={'initial','final090','final070'} and set(z.library)==set(LIBRARIES)
            clusters=pd.read_csv(score/'clusters.csv',dtype=str)
            assert clusters.cluster.nunique()==int(sm['n_clusters'])
            z['n_cells']=len(clusters);z['n_clusters']=int(sm['n_clusters'])
            z['n_hdbscan_noise']=int(clusters.cluster.eq('0').sum()) if cfg['method']=='HDBSCAN_R' else 0
            assert z[['macroF1_present','accuracy','coverage']].map(np.isfinite).all().all()
            neighbors.append(z);proof.append(dict(task='neighbors',path=str(p),n_rows=len(z)))
            artifacts[str(p/'manifest.json')]=sha(p/'manifest.json')
        else:
            base=OUT/'learning'/cfg['sample']/cfg['budget']/f"seed{cfg['model_seed']}"
            for library in LIBRARIES:
                p=base/library
                if not p.exists():
                    assert library==LIBRARIES[1];continue
                assert checked(p,'manifest.json','LEARNING_COMPLETE')
                m=json.loads((p/'manifest.json').read_text());tr=json.loads((p/'training_manifest.json').read_text())
                assert sha(p/'training_manifest.json')==m['training_manifest_sha256']
                for name,digest in m['files'].items():assert sha(p/name)==digest
                for name,digest in tr['outputs'].items():assert sha(p/name)==digest
                z=pd.read_csv(p/'metrics.csv');metrics.append(z)
                common=dict(sample=cfg['sample'],budget=cfg['budget'],model_seed=cfg['model_seed'],library=library,backup=library!=LIBRARIES[0],dl_status=tr['dl_status'],training_executed=tr['training_executed'],n_known=tr['n_known'],n_pool=tr['n_pool'],n_training_classes=tr['n_training_classes'],n_train=tr.get('n_train',0),n_validation=tr.get('n_validation',0))
                learning_states.append(common)
                if tr['training_executed']:
                    h=pd.read_csv(p/'learning_history.csv');assert list(h.epoch)==list(range(1,31))
                    assert (h.n==tr['n_train']).all() and (h.validation_n==tr['n_validation']).all()
                    assert np.isfinite(h.mean_loss).all() and h.accuracy.between(0,1).all()
                    if tr['n_validation']:
                        assert np.isfinite(h.validation_loss).all() and h.validation_accuracy.between(0,1).all()
                    assert len(z)==8 and set(z.epochs)=={5,10,20,30}
                    for cp in tr['checkpoints']:
                        d=p/f"epoch{cp['epoch']:02d}"
                        assert sha(d/'terminal.npz')==cp['terminal_sha256'] and sha(d/'model_state.pt')==cp['model_sha256']
                        pred=pd.read_csv(d/'predictions.csv.gz',dtype=str,keep_default_na=False)
                        cells=pd.read_csv(OLD/'GBM'/cfg['sample']/cfg['budget']/'cells.csv',dtype=str,keep_default_na=False)
                        assert list(pred.cell_id)==list(cells.cell_id)
                        with np.load(d/'terminal.npz',allow_pickle=False) as checkpoint:
                            for field in ['initial','final090','final070']:assert np.array_equal(pred[field].to_numpy(),checkpoint[field])
                        artifacts[str(d/'predictions.csv.gz')]=sha(d/'predictions.csv.gz')
                    for k,v in common.items():h[k]=v
                    histories.append(h)
                else:assert len(z)==2 and not (p/'training_history.json').exists()
                for name in ['clustering_and_terminal.png','clustering_and_terminal.pdf']:assert (p/name).stat().st_size>1000
                proof.append(dict(task='learning',path=str(p),n_rows=len(z)))
                artifacts[str(p/'manifest.json')]=sha(p/'manifest.json')
    dest=OUT/'summary';dest.mkdir(exist_ok=True)
    nt=pd.concat(neighbors,ignore_index=True);mt=pd.concat(metrics,ignore_index=True);ht=pd.concat(histories,ignore_index=True);st=pd.DataFrame(learning_states)
    nt.to_csv(dest/'neighbor_metrics.csv',index=False);mt.to_csv(dest/'checkpoint_metrics.csv',index=False);ht.to_csv(dest/'validation_histories.csv',index=False);st.to_csv(dest/'learning_status.csv',index=False)
    pd.DataFrame(proof).to_csv(dest/'verified_outputs.csv',index=False)
    variability=mt.groupby(['sample','budget','library','epochs','stage'],dropna=False).agg(n_seeds=('model_seed','nunique'),macroF1_mean=('macroF1_present','mean'),macroF1_sd=('macroF1_present','std'),macroF1_min=('macroF1_present','min'),macroF1_max=('macroF1_present','max'),coverage_mean=('coverage','mean'),coverage_min=('coverage','min'),coverage_max=('coverage','max')).reset_index()
    variability.to_csv(dest/'checkpoint_seed_variability.csv',index=False)
    anchors=nt[(nt.snn_k==20)&(nt.umap_neighbors==30)][['sample','budget','route','library','stage','macroF1_present','coverage']].rename(columns={'macroF1_present':'anchor_macroF1','coverage':'anchor_coverage'})
    delta=nt.merge(anchors,on=['sample','budget','route','library','stage'],how='left',validate='many_to_one')
    assert delta.anchor_macroF1.notna().all()
    delta['delta_macroF1']=delta.macroF1_present-delta.anchor_macroF1;delta['delta_coverage']=delta.coverage-delta.anchor_coverage
    delta.to_csv(dest/'neighbor_changes_from_anchor.csv',index=False)
    fig,axs=plt.subplots(3,2,figsize=(12,11),layout='constrained')
    for i,sample in enumerate(PILOTS):
        for j,budget in enumerate(BUDGETS):
            ax=axs[i,j]
            for seed,color in zip([0,1,42],['#0072B2','#D55E00','#009E73']):
                g=ht[(ht['sample']==sample)&ht.budget.eq(budget)&ht.model_seed.eq(seed)&ht.library.eq(LIBRARIES[0])]
                if len(g):
                    ax.plot(g.epoch,g.mean_loss,color=color,label=f'Seed {seed}: train')
                    ax.plot(g.epoch,g.validation_loss,color=color,ls='--',label=f'Seed {seed}: validation')
            if not len(ax.lines):ax.text(.5,.5,'No trainable primary refinement',transform=ax.transAxes,ha='center')
            classes=st[st['sample'].eq(sample)&st.budget.eq(budget)&st.library.eq(LIBRARIES[0])].n_training_classes.unique();assert len(classes)==1
            ax.axvline(10,color='#888',ls=':');ax.set(title=f'{sample} / {budget} / {classes[0]} known seed classes',xlabel='Epoch',ylabel='Softmax + cross-entropy loss')
    handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=3,fontsize=8)
    fig.suptitle('Real per-epoch training and held-out pseudo-label validation\nOriginal CM2 glioma context; validation does not supply independent biological truth')
    for ext in ['png','pdf']:fig.savefig(dest/f'validation_loss_summary.{ext}',dpi=160,bbox_inches='tight')
    plt.close(fig)
    # Every neighbour point is a terminal DL result; lines link prespecified parameter candidates only.
    fig,axs=plt.subplots(3,2,figsize=(12,11),layout='constrained')
    primary=nt[nt.library.eq(LIBRARIES[0])&nt.stage.eq('final090')]
    for i,sample in enumerate(PILOTS):
        for j,budget in enumerate(BUDGETS):
            ax=axs[i,j];g=primary[primary['sample'].eq(sample)&primary.budget.eq(budget)]
            for route,kind,key,color in [('PCA30_SNN','SNN k','snn_k','#0072B2'),('UMAP2_SNN','SNN k','snn_k','#009E73'),('UMAP2_SNN','UMAP neighbors','umap_neighbors','#D55E00'),('UMAP2_HDBSCAN_R','UMAP neighbors','umap_neighbors','#CC79A7')]:
                z=g[g.route.eq(route)&(g.umap_neighbors.eq(30) if key=='snn_k' else g.snn_k.eq(20))].sort_values(key)
                ax.plot(z[key],z.macroF1_present,'o-',color=color,label=route+' / '+kind)
            ax.set(title=f'{sample} / {budget}',xlabel='Prespecified neighbors',ylabel='All-cell terminal macro-F1',ylim=(0,1))
    handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=2,fontsize=8)
    fig.suptitle('Neighbor sensitivity: three size-selected GBM pilots\nFixed marker CM2_glioma_other / mean; fixed DL seed42, epoch10, threshold0.90')
    for ext in ['png','pdf']:fig.savefig(dest/f'neighbor_sensitivity_summary.{ext}',dpi=160,bbox_inches='tight')
    plt.close(fig)
    note=f'''# GBM reviewer B: neighbors and real validation curves

Completed:66 unique neighbor conditions on TKU4163, NL022, SN040 × HVG2000/5000. Two fixed marker contexts per partition. SNN k=10/20/40 in PCA30 and UMAP2; UMAP n.neighbors=15/30/60 for SNN and HDBSCAN. The SNN k20 / UMAP n30 duplicate is counted once. Default branches reuse verified original outputs after an independent three-route default-parity pilot. Every condition retains clustering and terminal DL figures.

Primary learning conditions:18 (three pilots × two budgets × seeds0/1/42), with {int(st.backup.sum())} separately identified backup-marker conditions. Real validation uses fixed split42 held-out marker pseudo-labels.30epochs are recorded; checkpoints5/10/20/30 retain models, probabilities and terminal predictions. All trainable10epoch checkpoints exactly match the existing original-width controls for their matching seeds. Original10epoch endpoints remain the reference. No biological truth was consulted during fitting, stopping or marker selection.

Scope: conditional three-pilot sensitivity, not an all-GBM optimum claim. Held-out pseudo-label agreement cannot independently validate cell biology. No-op conditions have no invented validation curve. Single-known-class training is separately identified; perfect pseudo-label accuracy in those fits does not establish multiclass discrimination. HDBSCANminPts50, SNNresolution0.5, UMAPcosine/min.dist0.3, MLP256/128, historicalSoftmax+CrossEntropy/Adamax and thresholds0.90/0.70 remain explicit. PCA and UMAP geometry use the original R objects and scorer; normalized selected RNA remains DL input. Fixed HVG2000 plot coordinates are display coordinates only.

PTC controls remain outside this GBM result package until the parent GBM gate is met.
'''
    (dest/'INTERPRETATION.md').write_text(note)
    write_json(dest/'manifest.json',dict(status='completed',n_neighbor_conditions=66,n_nominal_neighbor_grid_entries=72,n_primary_learning_conditions=18,n_backup_learning_conditions=int(st.backup.sum()),n_training_curves=int(st.training_executed.sum()),n_epoch_rows=len(ht),n_neighbor_metric_rows=len(nt),n_checkpoint_metric_rows=len(mt),
        original_source_unchanged=True,numerical_sources_match_parity=True,final_source_sha256={p.name:sha(p) for p in HERE.iterdir() if p.suffix in ['.py','.R']},parity_sha256=sha(OUT/'verification/parity.json'),verified_artifacts=artifacts,files={p.name:sha(p) for p in dest.iterdir() if p.is_file() and p.name not in ['manifest.json','COMPLETE']},job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest)
    state.update(status='completed',updated_at=utc(),completed=84,remaining=0,evidence=[str(dest/'manifest.json'),str(dest/'INTERPRETATION.md'),str(dest/'validation_loss_summary.png'),str(dest/'neighbor_sensitivity_summary.png')]);write_json(OUT/'status.json',state)
    print(json.dumps({k:v for k,v in json.loads((dest/'manifest.json').read_text()).items() if k not in ['verified_artifacts','files']}),flush=True)

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='refresh':print(json.dumps(refresh_status()))
    else:aggregate()
