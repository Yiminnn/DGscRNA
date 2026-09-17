import json,os,sys
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42})
    candidates=[OUT/'PTC_archived_CCA2000',*sorted((OUT/'PTC_ablation').iterdir())]
    ready=[p for p in candidates if (p/'evaluation/ABSTENTION_COMPLETE').exists()]
    if '--available' not in sys.argv:assert len(ready)==30,len(ready)
    aliases={'seurat_clusters':'P-SNN','seurat.UMAP_clusters':'U-SNN','hdbscan_clusters':'P-HDB','hdbscan.UMAP_clusters':'U-HDB',
       'PCA30_SNN':'P-SNN','PCA30_HDBSCAN_R':'P-HDB','UMAP2_SNN':'U-SNN','UMAP2_HDBSCAN_R':'U-HDB'}
    for prep in ready:
        f=pd.read_csv(prep/'evaluation/abstention_aware_metrics.csv.gz')
        f=f[f.stage.eq('final090')&f.endpoint.eq('strict_T_name_rule')]
        dest=prep/'figures';dest.mkdir(exist_ok=True)
        for group in ['NMT','TTU']:
            subset=f[f.scope.eq(group)]
            if subset.empty:continue
            piv=subset.pivot(index='library',columns=['route','cutoff'],values='macro_F1_T_nonT_unknown_as_error')
            fig,ax=plt.subplots(figsize=(7.2,max(3.4,.26*len(piv)+1.5)))
            im=ax.imshow(piv.to_numpy(),vmin=0,vmax=1,cmap='cividis',aspect='auto')
            ax.set_yticks(range(len(piv)));ax.set_yticklabels([x.replace('CellMarker_','CM ').replace('_',' ') for x in piv.index],fontsize=7)
            ax.set_xticks(range(len(piv.columns)));ax.set_xticklabels([aliases[a]+'\n'+b for a,b in piv.columns],fontsize=6)
            for i in range(len(piv)):
                for j in range(len(piv.columns)):
                    v=piv.iloc[i,j];ax.text(j,i,f'{v:.2f}',ha='center',va='center',fontsize=5.8,color='white' if v<.5 else 'black')
            fig.colorbar(im,ax=ax,fraction=.035,pad=.02,label='T / non-T macro F1; Unknown retained as error')
            ax.set_title(f'{prep.name} | {group}\nTerminal DL; strict T-name / TCR proxy; all marker contexts',fontsize=9)
            fig.tight_layout()
            for ext in ['png','pdf']:fig.savefig(dest/f'terminal_marker_grid_abstention_{group}.{ext}',dpi=300,bbox_inches='tight',facecolor='white')
            plt.close(fig)
        (dest/'ABSTENTION_PLOTTED').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'status':'complete'})+'\n')
        print('ABSTENTION_PLOTTED',prep.name,flush=True)

if __name__=='__main__':run()
