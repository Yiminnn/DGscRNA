"""Finite dimension/parameter and embedding-seed controls with full terminal outputs."""
import json
import os
from pathlib import Path
import subprocess
import sys
from common import OUT, RSCRIPT, L1, require_slurm, checked, complete, sha, write_json, utc
LIBRARIES=['CM2_glioma_other','CM2_primary_all_context']

def configurations(sample,budget,include_defaults=False):
    result=[]
    def add(space,method,minpts=50,res=.5,seed=42,kind='parameter'):
        name=f'{space}_{method}_minPts{minpts}_r{res}_seed{seed}'
        result.append(dict(name=name,sample=sample,budget=budget,space=space,method=method,
            minPts=minpts,resolution=res,embedding_seed=seed,kind=kind))
    if include_defaults:
        for space in ['PCA30','UMAP2']:
            for method in ['SNN','HDBSCAN_R']:add(space,method,kind='default_parity')
        return result
    for space in ['UMAP10','UMAP30','RNA_noDR']:
        for method in ['SNN','HDBSCAN_R']:add(space,method,kind='dimension')
    for space in ['PCA30','UMAP2']:
        for minpts in [25,100]:add(space,'HDBSCAN_R',minpts=minpts)
        for res in [.25,1.]:add(space,'SNN',res=res)
        for seed in [0,1,2,3]:
            for method in ['SNN','HDBSCAN_R']:add(space,method,seed=seed,kind='embedding_seed')
    assert len(result)==30
    return result

def run(cfg,parity=False):
    require_slurm()
    import numpy as np
    import pandas as pd
    import torch
    from sklearn.metrics import precision_recall_fscore_support
    import terminal
    torch.set_num_threads(4)
    sample=cfg['sample'];budget=cfg['budget'];name=cfg['name']
    dest=OUT/'GBM_representation_controls'/sample/budget/name;dest.mkdir(parents=True,exist_ok=True)
    if checked(dest):return
    source=Path(__file__).resolve().parent
    if not parity:
        proof=json.loads((OUT/'verification/representation_default_parity.json').read_text())
        assert proof['status']=='passed'
        assert proof['score_source_sha256']==sha(source/'score_representation_R.R')
    conf=dest/'config.json'
    if conf.exists():assert json.loads(conf.read_text())==cfg
    else:write_json(conf,cfg)
    subprocess.run([RSCRIPT,str(source/'prepare_representation_R.R'),str(conf)],check=True)
    subprocess.run([RSCRIPT,str(source/'score_representation_R.R'),sample,str(dest),str(conf)],check=True)
    route=cfg['space']+'_'+cfg['method'];score=dest/route
    sm=json.loads((score/'score_manifest.json').read_text())
    for aid,a in sm['arms'].items():
        if a['library'] in LIBRARIES and a['cutoff']=='mean':terminal.finish_route(score,only_arm=aid)
    if parity:
        original=OUT/'GBM'/sample/budget/route
        for file in ['clusters.csv','initial_calls.csv.gz']:
            a=pd.read_csv(score/file,dtype=str,keep_default_na=False)
            b=pd.read_csv(original/file,dtype=str,keep_default_na=False)
            assert a.equals(b),(name,file)
        for aid,a in sm['arms'].items():
            if a['library'] in LIBRARIES and a['cutoff']=='mean':
                got=np.load(score/'terminal'/aid/'terminal.npz',allow_pickle=False)
                ref=np.load(original/'terminal'/aid/'terminal.npz',allow_pickle=False)
                for key in ['initial','final090','final070','probabilities']:assert np.array_equal(got[key],ref[key]),(name,aid,key)
    # Evaluation follows frozen native predictions.
    truth=pd.read_csv(OUT/'evaluation_inputs'/sample/'truth.csv.gz',dtype=str,keep_default_na=False)
    mapping=pd.read_csv(OUT/'markers/panel_L1_mapping.csv',dtype=str,keep_default_na=False)
    maps={lib:dict(zip(g.panel,g.L1)) for lib,g in mapping.groupby('library')}
    y=truth.L1.to_numpy();rows=[];calls={}
    for aid,a in sm['arms'].items():
        if a['library'] not in LIBRARIES or a['cutoff']!='mean':continue
        td=score/'terminal'/aid;pred=pd.read_csv(td/'predictions.csv.gz',dtype=str,keep_default_na=False)
        assert np.array_equal(pred.cell_id,truth.cell_id)
        tm=json.loads((td/'terminal_manifest.json').read_text())
        for stage in ['initial','final090','final070']:
            p=np.asarray([maps[a['library']].get(v,'Unknown' if v in ['Unknown','Undecided',''] else 'UNMAPPABLE') for v in pred[stage]])
            _,_,f1,support=precision_recall_fscore_support(y,p,labels=L1,zero_division=0)
            rows.append(dict(**cfg,route=route,library=a['library'],stage=stage,macroF1_present=float(f1[support>0].mean()),
                accuracy=float((y==p).mean()),coverage=float((~pred[stage].isin(['Unknown','Undecided',''])).mean()),dl_status=tm['dl_status']))
            if stage=='final090':calls[a['library']]=p
    pd.DataFrame(rows).to_csv(dest/'metrics.csv',index=False)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    xy=pd.read_csv(OUT/'GBM'/sample/'hvg2000/UMAP2.csv',index_col=0)
    assert list(xy.index)==list(truth.cell_id)
    cl=pd.read_csv(score/'clusters.csv',dtype=str)
    colors={v:plt.get_cmap('tab20')(i) for i,v in enumerate(L1)}
    colors.update(Unknown='#bdbdbd',UNMAPPABLE='#7b614c',AMBIGUOUS_NEURON='#8c8c33',NO_L1_COUNTERPART='#7b614c')
    fig,axs=plt.subplots(2,2,figsize=(10,9),layout='constrained')
    panels=[('Author L1',y),('Clustering in '+cfg['space'],cl.cluster.to_numpy())]+list(calls.items())
    for i,(ax,(title,labs)) in enumerate(zip(axs.flat,panels)):
        unique=sorted(set(labs));palette={v:plt.get_cmap('tab20')(j%20) for j,v in enumerate(unique)} if i==1 else colors
        if i==1 and cfg['method']=='HDBSCAN_R':palette['0']='#bdbdbd'
        for lab in unique:
            keep=np.asarray(labs)==lab
            ax.scatter(xy.iloc[keep,0],xy.iloc[keep,1],s=2,c=[palette.get(lab,'#7b614c')],linewidths=0,rasterized=True)
        ax.set_title(title,fontsize=9);ax.set_xticks([]);ax.set_yticks([])
        ax.set(xlabel='Fixed HVG2000 display UMAP1',ylabel='Fixed display UMAP2')
    fig.suptitle(f'{sample} | {budget}\n{name}\nActual clustering space shown in title; all panels share original display coordinates',fontsize=11)
    fig.legend(handles=[Line2D([],[],marker='o',color='none',markerfacecolor=colors[v],label=v,markersize=4) for v in L1+['Unknown','UNMAPPABLE','AMBIGUOUS_NEURON']],
        loc='outside lower center',ncol=5,fontsize=7,frameon=False)
    for ext in ['png','pdf']:fig.savefig(dest/f'clustering_and_terminal.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)
    write_json(dest/'manifest.json',dict(status='completed',config=cfg,default_parity=parity,
        scope='Prespecified min/median/max-size samples. Conditional algorithm stability, not additional biological replicates or proof of all-cohort dimension optimality.',
        fixed_DL_and_full_RNA_scoring_within_budget=True,MLP_seed=42,split_seed=42,
        files={p.name:sha(p) for p in dest.iterdir() if p.suffix in ['.csv','.png','.pdf']},
        score_source_sha256=sha(source/'score_representation_R.R'),job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest)

if __name__=='__main__':
    if sys.argv[1]=='parity':
        cases=configurations('TKU4163','hvg2000',include_defaults=True)
        for cfg in cases:run(cfg,parity=True)
        write_json(OUT/'verification/representation_default_parity.json',dict(status='passed',cases=cases,
            exact=['partitions','all48marker seed conditions','two fixed-context terminal predictions and class probabilities'],
            score_source_sha256=sha(Path(__file__).resolve().parent/'score_representation_R.R'),job=os.environ['SLURM_JOB_ID']))
    elif sys.argv[1]=='tasklist':
        cfg=json.loads(Path(sys.argv[2]).read_text())[int(os.environ['SLURM_ARRAY_TASK_ID'])];run(cfg)
    else:run(json.loads(Path(sys.argv[1]).read_text()))
