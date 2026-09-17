"""Export every clustering branch and the complete terminal marker grid."""
import json,os,sys,textwrap
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'
from evaluation_rules import broad

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,
       'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42,'ps.fonttype':42})
    prep=Path(sys.argv[1]);dataset=sys.argv[2]
    dest=prep/'figures';dest.mkdir(exist_ok=True)
    assert (prep/'evaluation/COMPLETE').exists()
    routes=sorted(prep.glob('*/score_manifest.json'))
    assert len(routes)==4
    meta=[json.loads(p.read_text()) for p in routes]
    cells=pd.read_csv(routes[0].parent/'cells.csv',keep_default_na=False,dtype=str)
    ids=cells.cell_id.to_numpy()
    emb=pd.read_csv(prep/'UMAP2.csv',index_col=0).loc[ids].to_numpy()[:,:2]
    is_ptc=dataset=='PTC'
    if is_ptc:
        ref=pd.read_csv(ROOT/'results/hvg_ptc_20260916_v1/ptc_paper_baseline/paper_baseline_reference.csv.gz',keep_default_na=False).set_index('cell_id').loc[ids]
        sys.path.insert(0,str(ROOT/'handoff/ptc_recovery_20260916'))
        from label_rules import broad_lineage
        from functools import lru_cache
        broad_lineage=lru_cache(maxsize=65536)(broad_lineage)
        truth=ref.paper_final_native.map(broad_lineage).to_numpy()
        batch=ref['sample'].to_numpy()
        label_name=lambda x:broad_lineage(str(x))
        unit=prep.name
    else:
        unit=prep.parent.name
        ref=pd.read_csv(OUT/'inputs'/dataset/unit/'evaluation_only.csv.gz',keep_default_na=False).set_index('cell_id').loc[ids]
        truth=np.array([broad(s,dataset) for s in ref.truth])
        batch=ref.donor.to_numpy()
        pm=pd.read_csv(OUT/'markers/native_panel_metadata.csv',keep_default_na=False)
        names=dict(zip(pm.native,pm.cell_name))
        label_name=lambda x:broad(names.get(str(x),str(x)),dataset)
    palette=['#0072B2','#E69F00','#009E73','#CC79A7','#56B4E9','#D55E00','#F0E442','#222222']
    shapes=['o','s','^','D','v','P','X','>']
    route_titles={'PCA30_SNN':'PCA + SNN','PCA30_HDBSCAN_R':'PCA + HDBSCAN',
      'UMAP2_SNN':'UMAP + SNN','UMAP2_HDBSCAN_R':'UMAP + HDBSCAN',
      'seurat_clusters':'PCA + SNN','seurat.UMAP_clusters':'UMAP + SNN',
      'hdbscan_clusters':'PCA + HDBSCAN','hdbscan.UMAP_clusters':'UMAP + HDBSCAN'}
    def scatter(ax,labels,title,cluster=False):
        labels=np.array(labels,dtype=str)
        values=sorted(set(labels),key=lambda x:int(x) if cluster and x.lstrip('-').isdigit() else x)
        for j,value in enumerate(values):
            mask=labels==value
            color='#b4b4b4' if value in ['Unknown','Undecided'] else palette[j%len(palette)]
            ax.scatter(emb[mask,0],emb[mask,1],s=.5 if len(ids)>20000 else 1.2,
                color=color,alpha=.65,linewidths=0,marker=shapes[(j//len(palette))%len(shapes)],rasterized=True)
            if cluster:
                center=np.median(emb[mask],axis=0)
                ax.text(*center,value,fontsize=6,ha='center',va='center',
                    path_effects=[pe.withStroke(linewidth=1.6,foreground='white')])
        ax.set_title(title,pad=5);ax.set_xticks([]);ax.set_yticks([])
        ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
        for spine in ax.spines.values():spine.set_visible(False)
        return values
    def save(fig,name):
        fig.savefig(dest/(name+'.png'),dpi=300,bbox_inches='tight',facecolor='white')
        fig.savefig(dest/(name+'.pdf'),dpi=300,bbox_inches='tight',facecolor='white')
        plt.close(fig)
    fig,axs=plt.subplots(2,3,figsize=(7.2,5.6))
    for i,path in enumerate(routes):
        source=path.parent
        if (source/'clusters.csv').exists():cl=pd.read_csv(source/'clusters.csv',dtype=str).set_index('cell_id').loc[ids,'cluster']
        else:cl=pd.read_csv(source/'cells.csv',dtype=str).set_index('cell_id').loc[ids,'cluster']
        name=route_titles[source.name]
        scatter(axs.flat[i],cl,f'{chr(65+i)}  {name}\n{cl.nunique()} groups',cluster=True)
    truth_values=scatter(axs.flat[4],truth,'E  Archived annotation' if is_ptc else 'E  Curated common lineages')
    batch_values=scatter(axs.flat[5],batch,'F  Sample / donor')
    for values,x,title in [(truth_values,.27,'E  Lineages'),(batch_values,.79,'F  Sample / donor')]:
        handles=[Line2D([],[],color='#b4b4b4' if value in ['Unknown','Undecided'] else palette[j%8],
             marker=shapes[(j//8)%len(shapes)],linestyle='',markersize=3,
             label=textwrap.fill(value.replace('unmapped:',''),22)) for j,value in enumerate(values)]
        fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(x,-.02),ncol=2,
                   fontsize=6,frameon=False,title=title,title_fontsize=7)
    fig.suptitle(f'{unit} | all {len(ids):,} cells',fontsize=11,y=1.01)
    fig.subplots_adjust(hspace=.32,wspace=.22)
    fig.text(.5,-.005,'Same saved UMAP coordinates; numbers identify clusters (HDBSCAN 0 = noise).',ha='center',fontsize=7)
    save(fig,'all_cluster_branches')
    # Prespecified primary context, never chosen by looking at a performance maximum.
    final_sets=[];initial_sets=[];choice_notes=[]
    for path,m in zip(routes,meta):
        source=path.parent
        if is_ptc:
            before=np.empty(len(ids),dtype=object);after=before.copy()
            for group,lib,cut in [('NMT','CellMarker_Thyroid','none'),('TTU','Pubmed_34663816','mean')]:
                mask=ref.group.eq(group).to_numpy()
                if not mask.any():continue
                aid=next(k for k,a in m['arms'].items() if a['library']==lib and a['cutoff']==cut)
                z=np.load(source/'terminal'/aid/'terminal.npz',allow_pickle=False)
                sourceids=pd.read_csv(source/'cells.csv',dtype=str).cell_id
                index=pd.Index(sourceids).get_indexer(ids);assert (index>=0).all()
                before[mask]=z['initial'][index][mask];after[mask]=z['final090'][index][mask]
            choice_notes.append('NMT Thyroid/none; TTU Pubmed/mean')
        else:
            aid=next(k for k,a in m['arms'].items() if a['library']=='CM2_primary_normal' and a['cutoff']=='mean')
            z=np.load(source/'terminal'/aid/'terminal.npz',allow_pickle=False)
            sourceids=pd.read_csv(source/'cells.csv',dtype=str).cell_id
            index=pd.Index(sourceids).get_indexer(ids);assert (index>=0).all()
            before=z['initial'][index];after=z['final090'][index]
            choice_notes.append('CM2 primary normal / mean (fixed illustration)')
        initial_sets.append(np.array([label_name(s) for s in before]))
        final_sets.append(np.array([label_name(s) for s in after]))
    # One shared category map across panels is essential for interpreting refinement.
    vocabulary=sorted(set(truth)|set(np.concatenate(initial_sets))|set(np.concatenate(final_sets)))
    colors={v:('#b4b4b4' if v in ['Unknown','Undecided'] else palette[i%8]) for i,v in enumerate(vocabulary)}
    fig,axs=plt.subplots(2,4,figsize=(7.2,5.1))
    for row,sets in enumerate([initial_sets,final_sets]):
        for col,labels in enumerate(sets):
            ax=axs[row,col]
            for j,v in enumerate(vocabulary):
                mask=labels==v
                if not mask.any():continue
                ax.scatter(emb[mask,0],emb[mask,1],s=.55 if len(ids)>20000 else 1.1,color=colors[v],
                   marker=shapes[(j//8)%len(shapes)],alpha=.65,linewidths=0,rasterized=True)
            unknown=np.isin(labels,['Unknown','Undecided']).mean()
            route=route_titles[routes[col].parent.name]
            ax.set_title(f'{route}\n{"Marker only" if row==0 else "Terminal DL"}\nUnknown {unknown:.1%}',fontsize=7)
            ax.set_xticks([]);ax.set_yticks([])
            for spine in ax.spines.values():spine.set_visible(False)
    handles=[Line2D([],[],color=colors[v],marker=shapes[(j//8)%len(shapes)],linestyle='',markersize=3,
                    label=textwrap.fill(v.replace('unmapped:',''),28)) for j,v in enumerate(vocabulary)]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,-.01),ncol=3,fontsize=6,frameon=False)
    fig.suptitle(unit+' | '+choice_notes[0],fontsize=9,y=1.02)
    fig.subplots_adjust(hspace=.30,wspace=.10)
    save(fig,'marker_to_terminal_fixed_context')
    metric=pd.read_csv(prep/'evaluation/metrics.csv.gz')
    f=metric[metric.stage.eq('final090')]
    if is_ptc:
        views=[(g,f[f.scope.eq(g)&f.endpoint.eq('strict_T_name_rule')],'F1_T','Strict T-name / TCR agreement') for g in ['NMT','TTU'] if f.scope.eq(g).any()]
    else:views=[('all',f[f.endpoint.eq('common_lineage')],'macro_F1','Common-lineage macro F1')]
    for group,frame,value,title in views:
        piv=frame.pivot(index='library',columns=['route','cutoff'],values=value)
        fig,ax=plt.subplots(figsize=(7.2,max(3.4,.26*len(piv)+1.5)))
        im=ax.imshow(piv.to_numpy(),vmin=0,vmax=1,cmap='cividis',aspect='auto')
        ax.set_yticks(range(len(piv)));ax.set_yticklabels([x.replace('CellMarker_','CM ').replace('CM2_','CM2 ').replace('_',' ') for x in piv.index],fontsize=7)
        def short_route(x):
            return x.replace('seurat.UMAP_clusters','U-SNN').replace('seurat_clusters','P-SNN').replace('hdbscan.UMAP_clusters','U-HDB').replace('hdbscan_clusters','P-HDB').replace('PCA30_SNN','P-SNN').replace('PCA30_HDBSCAN_R','P-HDB').replace('UMAP2_SNN','U-SNN').replace('UMAP2_HDBSCAN_R','U-HDB')
        ax.set_xticks(range(len(piv.columns)));ax.set_xticklabels([short_route(a)+'\n'+b for a,b in piv.columns],fontsize=6)
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                v=piv.iloc[i,j]
                if np.isfinite(v):ax.text(j,i,f'{v:.2f}',ha='center',va='center',fontsize=5.8,color='white' if v<.5 else 'black')
        fig.colorbar(im,ax=ax,fraction=.035,pad=.02,label=title)
        ax.set_title(f'{unit} | {group}\nAll marker contexts, cutoffs and clustering branches; terminal DL',fontsize=9)
        fig.tight_layout();save(fig,'terminal_marker_grid_'+group)
    (dest/'FIGURES_COMPLETE').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'figure_version':2,'all_four_clusterings_plotted':True,'all_cells_plotted':len(ids),'formats':['PNG 300 dpi','PDF']})+'\n')
    print('PLOTTED',unit,len(ids),flush=True)

if __name__=='__main__':run()
