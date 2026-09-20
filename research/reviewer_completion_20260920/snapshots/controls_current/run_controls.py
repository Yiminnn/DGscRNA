"""Frozen GBM neighbour sensitivity and real pseudo-label validation curves."""
from pathlib import Path
import os, sys, json, subprocess, shutil
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OLD_CODE=ROOT/'handoff/paper_claim_validation_20260917'
sys.path.insert(0,str(OLD_CODE))
from common import OUT as OLD, RSCRIPT, L1, require_slurm, checked, complete, sha, write_json, utc
OUT=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/controls'
LIBRARIES=['CM2_glioma_other','CM2_primary_all_context']
PILOTS=['TKU4163','NL022','SN040'];BUDGETS=['hvg2000','hvg5000']

def tasks():
    result=[]
    for sample in PILOTS:
        for budget in BUDGETS:
            configs=[]
            for space in ['PCA30','UMAP2']:
                for k in [10,20,40]:configs.append((space,'SNN',k,30,'SNN_k'))
            for n in [15,30,60]:
                for method in ['SNN','HDBSCAN_R']:
                    if n==30 and method=='SNN':continue # exact duplicate SNN k20/UMAP n30
                    configs.append(('UMAP2',method,20,n,'UMAP_neighbors'))
            assert len(configs)==11
            for space,method,k,n,kind in configs:
                result.append(dict(task='neighbors',sample=sample,budget=budget,space=space,method=method,
                    snn_k=k,umap_neighbors=n,minPts=50,resolution=.5,embedding_seed=42,kind=kind,
                    name=f'{space}_{method}_k{k}_n{n}'))
            for seed in [0,1,42]:
                result.append(dict(task='learning',sample=sample,budget=budget,model_seed=seed,epochs=30))
    assert len(result)==84
    return result

def arm(source,library):
    sm=json.loads((source/'score_manifest.json').read_text())
    aids=[aid for aid,a in sm['arms'].items() if a['library']==library and a['cutoff']=='mean']
    assert len(aids)==1
    return aids[0],sm

def mapping_truth(sample):
    import pandas as pd
    truth=pd.read_csv(OLD/'evaluation_inputs'/sample/'truth.csv.gz',dtype=str,keep_default_na=False)
    m=pd.read_csv(OLD/'markers/panel_L1_mapping.csv',dtype=str,keep_default_na=False)
    return truth,{lib:dict(zip(g.panel,g.L1)) for lib,g in m.groupby('library')}

def evaluate(pred,truth,mapping,stage,context):
    import numpy as np
    from sklearn.metrics import precision_recall_fscore_support
    calls=pred[stage]
    p=np.asarray([mapping.get(v,'Unknown' if v in ['Unknown','Undecided',''] else 'UNMAPPABLE') for v in calls])
    y=truth.L1.to_numpy();_,_,f1,support=precision_recall_fscore_support(y,p,labels=L1,zero_division=0)
    return dict(**context,stage=stage,macroF1_present=float(f1[support>0].mean()),accuracy=float((y==p).mean()),
        coverage=float((~np.isin(calls,['Unknown','Undecided',''])).mean())),p

def plot_panels(sample,budget,title,panels,dest):
    import numpy as np,pandas as pd,matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    xy=pd.read_csv(OLD/'GBM'/sample/'hvg2000/UMAP2.csv',index_col=0)
    truth,_=mapping_truth(sample);assert list(xy.index)==list(truth.cell_id)
    colors={v:plt.get_cmap('tab20')(i) for i,v in enumerate(L1)};colors.update(Unknown='#aaa',UNMAPPABLE='#796547',AMBIGUOUS_NEURON='#9c9835')
    cols=2;rows=(len(panels)+1)//2
    fig,axs=plt.subplots(rows,cols,figsize=(10,4.5*rows),squeeze=False,layout='constrained')
    for ax,(label,labs,cluster) in zip(axs.flat,panels):
        labs=np.asarray(labs);palette={v:plt.get_cmap('tab20')(i%20) for i,v in enumerate(sorted(set(labs)))} if cluster else colors
        if cluster:palette['0']='#bdbdbd'
        for value in sorted(set(labs)):
            keep=labs==value;ax.scatter(xy.iloc[keep,0],xy.iloc[keep,1],s=2,c=[palette.get(value,'#796547')],lw=0,rasterized=True)
        ax.set(title=label,xlabel='Fixed HVG2000 display UMAP1',ylabel='Fixed display UMAP2',xticks=[],yticks=[])
    for ax in list(axs.flat)[len(panels):]:ax.set_visible(False)
    fig.suptitle(f'{sample} / {budget}\n{title}\nFixed display coordinates; actual clustering input recorded separately',fontsize=10)
    fig.legend(handles=[Line2D([],[],marker='o',color='none',markerfacecolor=colors[v],label=v,ms=4) for v in L1+['Unknown','UNMAPPABLE']],loc='outside lower center',ncol=5,fontsize=7,frameon=False)
    for ext in ['png','pdf']:fig.savefig(dest/f'clustering_and_terminal.{ext}',dpi=150,bbox_inches='tight')
    plt.close(fig)

def verify_equal(got,ref,fields):
    import numpy as np
    with np.load(got,allow_pickle=False) as a,np.load(ref,allow_pickle=False) as b:
        for key in fields:assert np.array_equal(a[key],b[key]),(str(got),key)

def run_neighbors(cfg,parity=False):
    require_slurm()
    import numpy as np,pandas as pd,torch,terminal
    torch.set_num_threads(4);terminal.OUT=OUT
    dest=OUT/'neighbors'/cfg['sample']/cfg['budget']/cfg['name'];dest.mkdir(parents=True,exist_ok=True)
    if checked(dest):return
    if not parity:assert checked(OUT/'verification','parity.json','PARITY_COMPLETE')
    cfg=dict(cfg,dest=str(dest));write_json(dest/'config.json',cfg)
    route=cfg['space']+'_'+cfg['method'];old=OLD/'GBM'/cfg['sample']/cfg['budget']/route
    score=dest/route;default=cfg['snn_k']==20 and cfg['umap_neighbors']==30
    if default and not parity:
        score.mkdir(exist_ok=True)
        for name in ['score_manifest.json','SCORE_COMPLETE','cells.csv','initial_calls.csv.gz','clusters.csv']:
            if not (score/name).exists():os.link(old/name,score/name)
        sm=json.loads((old/'score_manifest.json').read_text())
        for library in LIBRARIES:
            aid,_=arm(old,library);target=score/'terminal'/aid
            if not target.exists():shutil.copytree(old/'terminal'/aid,target,copy_function=os.link)
        reuse=dict(reused_exact_native_default=True,source=str(old),source_manifest_sha256=sha(old/'score_manifest.json'))
    else:
        subprocess.run([RSCRIPT,str(HERE/'prepare_neighbors.R'),str(dest/'config.json')],check=True)
        subprocess.run([RSCRIPT,str(HERE/'score_neighbors.R'),cfg['sample'],str(dest),str(dest/'config.json')],check=True)
        sm=json.loads((score/'score_manifest.json').read_text())
        for library in LIBRARIES:
            aid,_=arm(score,library);terminal.finish_route(score,only_arm=aid)
        reuse=dict(reused_exact_native_default=False)
    if parity:
        for name in ['clusters.csv','initial_calls.csv.gz']:
            assert pd.read_csv(score/name,dtype=str,keep_default_na=False).equals(pd.read_csv(old/name,dtype=str,keep_default_na=False)),(route,name)
        for library in LIBRARIES:
            aid,_=arm(score,library)
            verify_equal(score/'terminal'/aid/'terminal.npz',old/'terminal'/aid/'terminal.npz',['initial','final090','final070','probabilities','train_indices','validation_indices'])
    truth,maps=mapping_truth(cfg['sample']);metrics=[];panels=[('Author L1',truth.L1,False)]
    cl=pd.read_csv(score/'clusters.csv',dtype=str,keep_default_na=False);assert list(cl.cell_id)==list(truth.cell_id)
    panels.append(('Cluster: '+route,cl.cluster,True))
    for library in LIBRARIES:
        aid,_=arm(score,library);td=score/'terminal'/aid
        tr=json.loads((td/'training_manifest.json').read_text());pred=pd.read_csv(td/'predictions.csv.gz',dtype=str,keep_default_na=False)
        assert list(pred.cell_id)==list(truth.cell_id)
        for stage in ['initial','final090','final070']:
            row,calls=evaluate(pred,truth,maps[library],stage,dict(sample=cfg['sample'],budget=cfg['budget'],control=cfg['name'],route=route,library=library,
                snn_k=cfg['snn_k'],umap_neighbors=cfg['umap_neighbors'],dl_status=tr['dl_status'],training_executed=tr['training_executed']))
            metrics.append(row)
            if stage=='final090':panels.append((library+' / final DL 0.90',calls,False))
    pd.DataFrame(metrics).to_csv(dest/'metrics.csv',index=False)
    plot_panels(cfg['sample'],cfg['budget'],cfg['name'],panels,dest)
    write_json(dest/'manifest.json',dict(status='completed',config=cfg,default_parity=parity,**reuse,files={n:sha(dest/n) for n in ['metrics.csv','clustering_and_terminal.png','clustering_and_terminal.pdf']},job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest)

def train_learning(cfg,library,parity=False):
    require_slurm()
    import numpy as np,pandas as pd,torch
    sys.path.insert(0,str(HERE));import refine_learning as refine
    torch.set_num_threads(4)
    sample=cfg['sample'];budget=cfg['budget'];seed=cfg['model_seed'];epochs=cfg['epochs']
    dest=OUT/('learning_parity' if parity else 'learning')/sample/budget/f'seed{seed}'/library;dest.mkdir(parents=True,exist_ok=True)
    if checked(dest,'manifest.json','LEARNING_COMPLETE'):return json.loads((dest/'training_manifest.json').read_text())
    source=OLD/'GBM'/sample/budget/'UMAP2_HDBSCAN_R';aid,sm=arm(source,library)
    assert checked(source,'score_manifest.json','SCORE_COMPLETE')
    cells=pd.read_csv(source/'cells.csv',dtype=str,keep_default_na=False)
    initial=pd.read_csv(source/'initial_calls.csv.gz',dtype=str,keep_default_na=False)
    assert list(cells.cell_id)==list(initial.cell_id)
    binary=Path(sm['DL_binary']);assert sha(binary)==sm['DL_binary_sha256']
    x=np.memmap(binary,mode='r',dtype='<f4',shape=(len(cells),int(sm['DL_features'])))
    refine.PARAMS.update(epochs=epochs,model_seed=seed,split_seed=42,input='Native R normalized selected RNA; original fixed genes/order')
    if not checked(dest,'training_manifest.json','COMPLETE'):
        refine.train_cache(x,initial[aid].to_numpy(dtype=str),dest,dict(source_score=str(source),score_sha256=sha(source/'score_manifest.json'),
            DL_sha256=sha(binary),library=library,cutoff='mean',seed_column=aid,reference_labels_used_for_fit=False,split_seed=42,source_sha256=sha(refine.__file__)))
    tr=json.loads((dest/'training_manifest.json').read_text())
    if parity:
        verify_equal(dest/'terminal.npz',source/'terminal'/aid/'terminal.npz',['initial','final090','final070','probabilities','train_indices','validation_indices'])
        if tr['training_executed']:
            oldhistory=json.loads((source/'terminal'/aid/'training_history.json').read_text());newhistory=json.loads((dest/'training_history.json').read_text())
            for a,b in zip(oldhistory,newhistory):
                for field in a:assert a[field]==b[field],field
    truth,maps=mapping_truth(sample);assert list(truth.cell_id)==list(cells.cell_id)
    metrics=[];panels=[]
    if tr['training_executed']:
        history=json.loads((dest/'training_history.json').read_text());pd.DataFrame(history).to_csv(dest/'learning_history.csv',index=False)
        for epoch in [5,10,20,30]:
            if epoch>epochs:continue
            td=dest/f'epoch{epoch:02d}';z=np.load(td/'terminal.npz',allow_pickle=False)
            pred=cells[['cell_id']].copy()
            for stage in ['initial','final090','final070']:pred[stage]=z[stage]
            pred.to_csv(td/'predictions.csv.gz',index=False)
            for stage in ['final090','final070']:
                row,calls=evaluate(z,truth,maps[library],stage,dict(sample=sample,budget=budget,library=library,model_seed=seed,epochs=epoch,dl_status=tr['dl_status'],training_executed=True))
                metrics.append(row)
                if stage=='final090':panels.append((f'Epoch {epoch} terminal DL 0.90',calls,False))
        # Independent saved endpoint parity at original epoch, all three seeds.
        name='original' if seed==42 else f'model_seed{seed}'
        ref=OLD/'GBM_DL_controls'/sample/budget/name/'UMAP2_HDBSCAN_R/terminal'/aid/'terminal.npz'
        verify_equal(dest/'epoch10/terminal.npz',ref,['initial','final090','final070','probabilities','train_indices','validation_indices'])
        import matplotlib
        matplotlib.use('Agg');import matplotlib.pyplot as plt
        fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
        h=pd.DataFrame(history)
        for ax,key,ylab in zip(axs,['loss','accuracy'],['Loss','Accuracy']):
            ax.plot(h.epoch,h['mean_loss' if key=='loss' else 'accuracy'],label='Train marker pseudo-labels')
            ax.plot(h.epoch,h['validation_'+key],label='Held-out marker pseudo-labels')
            ax.axvline(10,color='#888',ls=':');ax.set(xlabel='Epoch',ylabel=ylab);ax.legend(fontsize=7)
        fig.suptitle(f'{sample} {budget} {library}, model seed {seed}, split42\n{tr["n_training_classes"]} known seed classes; pseudo-label validation is not biological validation',fontsize=10)
        for ext in ['png','pdf']:fig.savefig(dest/f'learning_curves.{ext}',dpi=160,bbox_inches='tight')
        plt.close(fig)
    else:
        z=np.load(dest/'terminal.npz',allow_pickle=False)
        for stage in ['final090','final070']:
            row,calls=evaluate(z,truth,maps[library],stage,dict(sample=sample,budget=budget,library=library,model_seed=seed,epochs=0,dl_status=tr['dl_status'],training_executed=False));metrics.append(row)
            if stage=='final090':panels=[('No training: '+tr['dl_status'],calls,False)]
    pd.DataFrame(metrics).to_csv(dest/'metrics.csv',index=False)
    plot_panels(sample,budget,f'{library}: checkpoints / seed{seed}',panels,dest)
    write_json(dest/'manifest.json',dict(status='completed',config=cfg,library=library,backup=library!=LIBRARIES[0],training_executed=tr['training_executed'],dl_status=tr['dl_status'],
        validation_target='Held-out marker pseudo-labels, not independent biological truth',original10epoch_parity=True,
        training_manifest_sha256=sha(dest/'training_manifest.json'),files={p.name:sha(p) for p in dest.iterdir() if p.suffix in ['.csv','.png','.pdf']},job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest,'manifest.json','LEARNING_COMPLETE');return tr

def run_learning(cfg):
    assert checked(OUT/'verification','parity.json','PARITY_COMPLETE')
    tr=train_learning(cfg,LIBRARIES[0])
    if not tr['training_executed']:train_learning(cfg,LIBRARIES[1])

def parity():
    configs=[c for c in tasks() if c['task']=='neighbors' and c['sample']=='TKU4163' and c['budget']=='hvg2000' and c['snn_k']==20 and c['umap_neighbors']==30]
    assert len(configs)==3
    for cfg in configs:run_neighbors(cfg,parity=True)
    train_learning(dict(task='learning',sample='TKU4163',budget='hvg2000',model_seed=42,epochs=10),LIBRARIES[0],parity=True)
    p=OUT/'verification';p.mkdir(exist_ok=True)
    write_json(p/'parity.json',dict(status='passed',partitions=3,marker_seed_arms_per_partition=48,terminal_contexts_per_partition=2,
        exact=['clusters','all marker seed labels','terminal labels','probabilities','training/validation cell split','10-epoch training histories'],
        source_hashes={p.name:sha(p) for p in HERE.iterdir() if p.suffix in ['.py','.R']},job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(p,'parity.json','PARITY_COMPLETE')

if __name__=='__main__':
    if sys.argv[1]=='freeze':
        OUT.mkdir(parents=True,exist_ok=True);write_json(OUT/'tasks.json',tasks())
        write_json(OUT/'status.json',dict(stage='B_GBM',status='preparing_parity',updated_at=utc(),jobs=[],completed=0,remaining=84,evidence=[str(OUT/'tasks.json')]))
    elif sys.argv[1]=='parity':parity()
    elif sys.argv[1]=='task':
        cfg=json.loads((OUT/'tasks.json').read_text())[int(os.environ['SLURM_ARRAY_TASK_ID'])]
        try:
            (run_neighbors if cfg['task']=='neighbors' else run_learning)(cfg)
        except Exception as e:
            import traceback
            write_json(OUT/'failures'/f'task_{os.environ["SLURM_ARRAY_TASK_ID"]}.json',dict(config=cfg,error=str(e),traceback=traceback.format_exc(),job=os.environ['SLURM_JOB_ID'],at=utc()))
            raise
        finally:
            from control_summary import refresh_status
            refresh_status()
    elif sys.argv[1]=='aggregate':
        from control_summary import aggregate
        aggregate()
