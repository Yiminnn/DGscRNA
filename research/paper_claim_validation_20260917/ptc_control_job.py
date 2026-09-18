"""Separate PTC preparation, scoring, refinement and evaluation SLURM stages."""
import json
import os
from pathlib import Path
import subprocess
import sys
from common import OLD,RSCRIPT,sha,checked,write_json,complete,utc
from ptc_followup_common import PTC,ANCHORS,CONTROL_ROUTES,require_ptc,read_reference,original_arm

SOURCE=Path(__file__).resolve().parent

def context_arms(cfg,route):
    source=Path(cfg['dest'])/route
    if cfg['kind']=='retention':return list(json.loads((source/'score_manifest.json').read_text())['arms'])
    if cfg['kind']=='default_parity':return [original_arm(source,cfg['group'])]
    contexts=json.loads((PTC/'selection/frozen_seed_contexts.json').read_text())
    old=OLD/'PTC_ablation'/f"PTC_{cfg['group']}_CCA{cfg['budget']}"/route
    return contexts[str(old)]['arm_ids']

def score(cfg,route):
    require_ptc()
    import numpy as np
    import pandas as pd
    dest=Path(cfg['dest']);source=dest/route
    assert checked(dest,'prepare_manifest.json','PREPARED')
    env=dict(os.environ,DGSCRNA_ONLY_ROUTE=route)
    subprocess.run([RSCRIPT,str(SOURCE/'ptc_score_R.R'),cfg['name'],str(dest)],env=env,check=True)
    assert checked(source,'score_manifest.json','SCORE_COMPLETE')
    if cfg['kind']=='default_parity':
        original=OLD/'PTC_ablation'/f"PTC_{cfg['group']}_CCA{cfg['budget']}"/route
        for name in ['clusters.csv','initial_calls.csv.gz']:
            left=pd.read_csv(source/name,dtype=str,keep_default_na=False)
            right=pd.read_csv(original/name,dtype=str,keep_default_na=False)
            assert left.equals(right),(str(original),name)
        write_json(source/'default_parity.json',dict(status='passed',reference=str(original),
            exact=['all cell partitions','all17libraries x3cutoffs initial calls'],
            score_source_sha256=sha(SOURCE/'ptc_score_R.R'),job=os.environ['SLURM_JOB_ID']))
    if cfg['kind']=='retention':
        original=OLD/'PTC_ablation'/f"PTC_{cfg['group']}_GEOMETRY{cfg['budget']}_FIXED_CCAall_DL2000"/route
        left=pd.read_csv(source/'clusters.csv',dtype=str);right=pd.read_csv(original/'clusters.csv',dtype=str)
        assert left.equals(right)
        pm=json.loads((dest/'prepare_manifest.json').read_text())
        om=json.loads((original/'score_manifest.json').read_text())
        assert pm['DL_binary_sha256']==om['DL_binary_sha256']
        config=dest/'retention_audit_config.json'
        if config.exists():assert json.loads(config.read_text())==cfg
        else:write_json(config,cfg)
        subprocess.run([RSCRIPT,str(SOURCE/'ptc_retention_audit.R'),str(config),route],check=True)
        write_json(source/'retention_invariants.json',dict(status='passed',partition_cells_exact=True,
            same_DL_binary=True,reference=str(original),geometry_unchanged=True,
            shared_gene_DEG_audit_sha256=sha(source/'retention_DEG_audit.json'),job=os.environ['SLURM_JOB_ID']))

def terminal(cfg,route,aid):
    require_ptc()
    from terminal import threads
    from ptc_terminal import fit
    threads();source=Path(cfg['dest'])/route
    assert aid in context_arms(cfg,route)
    fit(source,aid,source/'terminal'/aid,42)

def mlp(cfg):
    require_ptc()
    from terminal import threads
    from ptc_terminal import fit
    from ptc_evaluate import evaluate_source
    import pandas as pd
    threads();dest=Path(cfg['dest']);dest.mkdir(parents=True,exist_ok=True)
    if checked(dest):return
    source=Path(cfg['source'])
    for aid in cfg['arm_ids']:
        fit(source,aid,dest/'terminal'/aid,cfg['seed'],
            parity_reference=source/'terminal'/aid if cfg['seed']==42 else None)
    rows,states=evaluate_source(source,read_reference(),{aid:dest/'terminal'/aid for aid in cfg['arm_ids']},
        context=dict(family='MLP_seed',model_seed=cfg['seed'],control_name=cfg['name']))
    frame=pd.DataFrame(rows);statuses=pd.DataFrame(states)
    frame=frame[frame.group==cfg['group']];statuses=statuses[statuses.group==cfg['group']]
    frame.to_csv(dest/'metrics.csv.gz',index=False);statuses.to_csv(dest/'terminal_statuses.csv.gz',index=False)
    write_json(dest/'manifest.json',dict(status='complete',config=cfg,representation_unchanged=True,
        default_parity=cfg['seed']==42,job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest)

def plot_control(cfg,route,ref):
    require_ptc()
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from ptc_label_rules import broad_lineage
    source=Path(cfg['dest'])/route;aid=original_arm(source,cfg['group'])
    pred=pd.read_csv(source/'terminal'/aid/'predictions.csv.gz',keep_default_na=False)
    metadata=ref.loc[pred.cell_id]
    xy=pd.read_csv(OLD/'PTC_ablation'/f"PTC_{cfg['group']}_CCA2000"/'UMAP2.csv',index_col=0).loc[pred.cell_id]
    cl=pd.read_csv(source/'clusters.csv',dtype=str).set_index('cell_id').loc[pred.cell_id,'cluster']
    mapped=lambda values:np.asarray([broad_lineage(str(v)) for v in values])
    panels=[('Saved paper lineage (concordance only)',mapped(metadata.paper_native)),
        ('Productive high-confidence TCR',np.where(metadata.TCR_cell_high_confidence_productive_TCR,'Detected','Not detected')),
        ('Clustering: '+route,cl.to_numpy()),('Fixed historical marker: initial',mapped(pred.initial)),
        ('Fixed historical marker: terminal DL',mapped(pred.final090))]
    lineage=['T','NK','NKT','B','Plasma','Myeloid','Endothelial','Stromal','Epithelial','Lymphoid_ambiguous','Tumor_unspecified','Other','Unknown']
    palette={v:plt.get_cmap('tab20')(i) for i,v in enumerate(lineage)}
    palette.update(Unknown='#bdbdbd',Detected='#d95f02',**{'Not detected':'#bdbdbd'})
    fig,axes=plt.subplots(2,3,figsize=(15,9),layout='constrained')
    for i,(ax,(title,labs)) in enumerate(zip(axes.flat,panels)):
        colors={v:plt.get_cmap('tab20')(j%20) for j,v in enumerate(sorted(set(labs)))} if i==2 else palette
        if i==2 and route.endswith('HDBSCAN_R'):colors['0']='#bdbdbd'
        for label in sorted(set(labs)):
            mask=labs==label
            ax.scatter(xy.iloc[mask,0],xy.iloc[mask,1],s=.5,c=[colors.get(label,'#777777')],linewidths=0,rasterized=True)
        ax.set_title(title,fontsize=10);ax.set_xticks([]);ax.set_yticks([])
        ax.set(xlabel='Fixed group CCA2000 display UMAP1',ylabel='Fixed display UMAP2')
    axes.flat[5].axis('off')
    axes.flat[5].legend(handles=[Line2D([],[],marker='o',color='none',markerfacecolor=palette[v],label=v,markersize=5)
        for v in lineage+['Detected','Not detected']],loc='center',ncol=2,fontsize=8,frameon=False)
    fig.suptitle(f"{cfg['name']}\n{route}; MLP seed42; retained initial calls; TCR absence does not establish non-T identity",fontsize=12)
    for ext in ['png','pdf']:fig.savefig(source/f'clustering_and_terminal.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)

def finish(cfg):
    require_ptc()
    import pandas as pd
    from ptc_evaluate import evaluate_source
    dest=Path(cfg['dest']);evaluation=dest/'evaluation';evaluation.mkdir(exist_ok=True)
    if checked(evaluation):return
    ref=read_reference();rows=[];statuses=[]
    for route in CONTROL_ROUTES:
        source=dest/route;arms=context_arms(cfg,route)
        assert all(checked(source/'terminal'/aid,'terminal_manifest.json','TERMINAL_COMPLETE') for aid in arms)
        r,s=evaluate_source(source,ref,{aid:source/'terminal'/aid for aid in arms},
            context=dict(family=cfg['kind'],representation_seed=cfg['seed'],model_seed=42,geometry_budget=cfg['budget']))
        rows.extend(r);statuses.extend(s);plot_control(cfg,route,ref)
    pd.DataFrame(rows).to_csv(evaluation/'metrics.csv.gz',index=False)
    pd.DataFrame(statuses).to_csv(evaluation/'terminal_statuses.csv.gz',index=False)
    write_json(evaluation/'manifest.json',dict(status='complete',config=cfg,
        marker_scope='all17x3' if cfg['kind']=='retention' else 'frozen original and training-patient-selected contexts',
        scoring_completed_for_all17x3=True,unrequested_terminal_arms_are_not_claimed_as_results=True,
        job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(evaluation)

if __name__=='__main__':
    require_ptc();mode=sys.argv[1];path=Path(sys.argv[2]);cfg=json.loads(path.read_text())
    if mode=='prepare':subprocess.run([RSCRIPT,str(SOURCE/'ptc_prepare_followup.R'),str(path)],check=True)
    elif mode=='score':score(cfg,sys.argv[3])
    elif mode=='terminal':terminal(cfg,sys.argv[3],sys.argv[4])
    elif mode=='MLP':mlp(cfg)
    elif mode=='finish':finish(cfg)
    else:raise ValueError(mode)
