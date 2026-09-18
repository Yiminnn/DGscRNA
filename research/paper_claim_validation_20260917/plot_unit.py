"""Every clustering route and final annotation on identical display coordinates."""
import json
import os
from pathlib import Path
import shutil
import sys
from common import OUT, L1, ROUTES, require_slurm, write_json, sha, complete, checked, utc

def run(prep, refresh_marker_legend=False):
    require_slurm()
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    prep=Path(prep);pm=json.loads((prep/'prepare_manifest.json').read_text())
    sample=pm['sample'];dest=prep/'figures';dest.mkdir(exist_ok=True)
    if checked(dest,'manifest.json','FIGURES_COMPLETE'):
        prior=json.loads((dest/'manifest.json').read_text())
        if not refresh_marker_legend or prior.get('all_marker_legends_included'):return
        backup=prep/'figures_before_marker_legend'
        if not backup.exists():shutil.copytree(dest,backup)
    reference=OUT/'GBM'/sample/'hvg2000'
    assert checked(reference,'prepare_manifest.json','PREPARED'),'Common display hvg2000 must finish first'
    assert checked(prep/'evaluation')
    t=pd.read_csv(OUT/'evaluation_inputs'/sample/'truth.csv.gz',dtype=str,keep_default_na=False)
    coords=pd.read_csv(reference/'UMAP2.csv',index_col=0)
    assert list(coords.index)==list(t.cell_id)
    coords=coords.to_numpy()
    cm=plt.get_cmap('tab20')
    colors={v:cm(i) for i,v in enumerate(L1)}
    colors.update(Unknown='#bdbdbd',UNMAPPABLE='#7b614c',AMBIGUOUS_NEURON='#8c8c33',NO_L1_COUNTERPART='#7b614c')
    mapping=pd.read_csv(OUT/'markers/panel_L1_mapping.csv',dtype=str,keep_default_na=False)
    maps={lib:dict(zip(g.panel,g.L1)) for lib,g in mapping.groupby('library')}
    def scatter(ax,labels,title,clustering=False,hdb=False):
        labs=np.asarray(labels,dtype=str)
        unique=sorted(set(labs),key=lambda v:int(v) if clustering else v)
        palette={v:cm(i%20) for i,v in enumerate(unique)} if clustering else colors
        if clustering and hdb:palette['0']='#bdbdbd'
        for v in unique:
            take=labs==v
            ax.scatter(coords[take,0],coords[take,1],s=2.2,c=[palette.get(v,'#7b614c')],
                       linewidths=0,rasterized=True,alpha=.8)
            if clustering and take.any():
                xy=np.median(coords[take],axis=0)
                ax.text(*xy,'noise' if hdb and v=='0' else v,fontsize=6,ha='center',bbox=dict(facecolor='white',alpha=.7,edgecolor='none',pad=.5))
        ax.set_title(title,fontsize=9);ax.set_xticks([]);ax.set_yticks([])
        ax.set_xlabel('Fixed HVG2000 UMAP1',fontsize=7);ax.set_ylabel('UMAP2',fontsize=7)
        for s in ax.spines.values():s.set_visible(False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(3,3,figsize=(12,10),layout='constrained')
    scatter(axes.flat[0],t.L1,'Original author L1 labels')
    metric=pd.read_csv(prep/'evaluation/metrics.csv',dtype={'cutoff':str})
    for j,route in enumerate(ROUTES):
        cl=pd.read_csv(prep/route/'clusters.csv',dtype=str)
        assert list(cl.cell_id)==list(t.cell_id)
        hdb='HDBSCAN' in route
        n_clusters=cl.loc[~cl.cluster.eq('0') if hdb else cl.cluster.notna(),'cluster'].nunique()
        description=f'{n_clusters} clusters'+(f'; noise {cl.cluster.eq("0").mean():.1%}' if hdb else '')
        scatter(axes.flat[j+1],cl.cluster,f'{route}\n{description}',True,hdb)
        pp=pd.read_csv(prep/route/'terminal/L00_mean/predictions.csv.gz',dtype=str,keep_default_na=False)
        assert list(pp.cell_id)==list(t.cell_id)
        labels=[maps['CM2_glioma_other'].get(v,'Unknown' if v in ['Unknown','Undecided'] else 'UNMAPPABLE') for v in pp.final090]
        m=metric[(metric.route==route)&(metric.library=='CM2_glioma_other')&(metric.cutoff=='mean')&(metric.stage=='terminal090')&(metric.family=='native_R_budget')].iloc[0]
        scatter(axes.flat[j+5],labels,f'{route}: final DL\nmacro-F1={m.macroF1_present:.3f}, coverage={m.coverage:.1%}')
    fig.suptitle(f'{sample} | {pm["budget"]} | Native R workflow\nFixed CM2_glioma_other / mean cutoff; display coordinates shared across all conditions',fontsize=13)
    handles=[Line2D([],[],marker='o',color='none',markerfacecolor=colors[v],label=v,markersize=5) for v in L1+['Unknown','UNMAPPABLE','AMBIGUOUS_NEURON']]
    fig.legend(handles=handles,loc='outside lower center',ncol=5,fontsize=7,frameon=False)
    for suffix in ['png','pdf']:fig.savefig(dest/f'all_cluster_routes.{suffix}',dpi=180,bbox_inches='tight')
    plt.close(fig)
    pilot=json.loads((OUT/'protocol/input_audit.json').read_text())
    if sample==pilot['display_sample']:
        libraries=json.loads((OUT/'markers/manifest.json').read_text())['libraries']
        for route in ROUTES:
            fig,axes=plt.subplots(len(libraries),2,figsize=(8,2.3*len(libraries)),layout='constrained')
            sm=json.loads((prep/route/'score_manifest.json').read_text())
            for i,lib in enumerate(libraries):
                aid=next(k for k,a in sm['arms'].items() if a['library']==lib and a['cutoff']=='mean')
                pp=pd.read_csv(prep/route/'terminal'/aid/'predictions.csv.gz',dtype=str,keep_default_na=False)
                for j,key in enumerate(['initial','final090']):
                    mapped=[maps[lib].get(v,'Unknown' if v in ['Unknown','Undecided'] else 'UNMAPPABLE') for v in pp[key]]
                    scatter(axes[i,j],mapped,lib+(' | marker-only ablation' if j==0 else ' | terminal DL'))
            fig.suptitle(f'{sample} | {pm["budget"]} | {route}\nAll marker contexts; fixed mean cutoff and display coordinates',fontsize=13)
            fig.legend(handles=handles,loc='outside lower center',ncol=3,fontsize=8,frameon=False)
            for suffix in ['png','pdf']:fig.savefig(dest/f'{route}_all_markers.{suffix}',dpi=130,bbox_inches='tight')
            plt.close(fig)
    files={p.name:sha(p) for p in dest.iterdir() if p.suffix in ['.png','.pdf']}
    write_json(dest/'manifest.json',dict(status='completed',sample=sample,budget=pm['budget'],clustering_routes=ROUTES,
        all_cells_plotted=len(t),display_coordinates=str(reference/'UMAP2.csv'),display_sha256=sha(reference/'UMAP2.csv'),
        all_marker_legends_included=sample==pilot['display_sample'],
        display_not_used_as_substitute_for_clustering_input=True,files=files,source_sha256=sha(__file__),
        job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest,'manifest.json','FIGURES_COMPLETE');print('FIGURES_COMPLETE',sample,pm['budget'],flush=True)

if __name__=='__main__':run(sys.argv[1])
