"""Patient/tissue composition, marker expression and Unknown/TCR-stratified QC."""
import json
import os
from pathlib import Path
from ptc_common import BASE,GROUPS,require_slurm,sha,utc,write_json
from label_rules import broad_lineage

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    dest=BASE/'biology_detail';dest.mkdir(exist_ok=True)
    figs=BASE/'figures';figs.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':170})
    ref=pd.read_csv(BASE/'evaluation_reference/reference_cells.csv.gz',index_col=0)
    genes=['CD3D','CD3E','TRAC','NKG7','GNLY','KLRD1','MS4A1','CD79A','MZB1','LYZ','LST1','EPCAM','KRT19','TG','COL1A1','CD4','CD8A','FOXP3']
    lineage_order=['T','NKT','NK','B','Plasma','Myeloid','Endothelial','Stromal','Epithelial','Lymphoid_ambiguous','Tumor_unspecified','Other','Unknown']
    cmap=plt.get_cmap('tab20');palette={l:cmap(i) for i,l in enumerate(lineage_order)};palette['Unknown']='#b9b9b9'
    qc=[];expression=[];composition=[];datasets={};inputs={};figure_rows=[]
    for group,samples in GROUPS.items():
      x=pd.read_csv(BASE/'biology'/f'{group}.gene_expression.csv.gz',index_col=0)
      assert set(genes)<=set(x.columns)
      for correction in ['NONE','CCA','HARMONY']:
        part=BASE/'pooled'/group/correction/'UMAP2_HDBSCAN_R';score=part/'score_RNA'
        d=score/'refinement/L07_mean'
        assert (d/'TERMINAL_COMPLETE').exists()
        tm=json.loads((d/'terminal_manifest.json').read_text())
        assert tm['terminal_sha256']==sha(d/'terminal.npz')
        cells=pd.read_csv(score/'initial_calls.csv.gz',usecols=['cell_id']).cell_id.to_list()
        r=ref.loc[cells].copy();xx=x.loc[cells,genes]
        z=np.load(d/'terminal.npz',allow_pickle=False)
        umap=pd.read_csv(part.parent/'UMAP2.csv',index_col=0).loc[cells]
        assert len(r)==len(z['final090']) and umap.shape==(len(r),2)
        inputs[str(d.relative_to(BASE))]=dict(terminal_sha256=sha(d/'terminal.npz'),cell_order_sha256=sha(score/'initial_calls.csv.gz'),UMAP_sha256=sha(part.parent/'UMAP2.csv'))
        for stage,key in [('marker_only_ablation','initial'),('terminal_DL090','final090')]:
          native=z[key];lookup={v:broad_lineage(v) for v in np.unique(native)}
          broad=np.asarray([lookup[v] for v in native]);assert set(broad)<=set(lineage_order)
          r['prediction']=broad
          r['prediction_category']=np.where(broad=='T','Called T',np.where(broad=='Unknown','Unknown','Other called'))
          detected=r.TCR_cell_high_confidence_productive_TCR.astype(bool)
          r['TCR_stratum']=np.where(detected,'TCR+','TCR not detected')
          r['stratum']=r.prediction_category+' / '+r.TCR_stratum
          for sample in samples:
            sub=r[r['sample'].eq(sample)];patient=sub.patient.iloc[0]
            values=sub.prediction.value_counts()
            for lineage in lineage_order:composition.append(dict(group=group,correction=correction,stage=stage,sample=sample,patient=patient,broad_lineage=lineage,n=int(values.get(lineage,0)),fraction=float(values.get(lineage,0)/len(sub))))
            for stratum,rr in sub.groupby('stratum'):
              qc.append(dict(group=group,correction=correction,stage=stage,sample=sample,patient=patient,stratum=stratum,n_cells=len(rr),
                fraction_of_sample=len(rr)/len(sub),median_nCount=rr.nCount_RNA.median(),median_nFeature=rr.nFeature_RNA.median(),
                median_percent_mt=rr['percent.mt'].median(),T_RNA_support_fraction=rr.T_RNA_support_ge2.mean()))
              for gene in genes:
                values=xx.loc[rr.index,gene]
                expression.append(dict(group=group,correction=correction,stage=stage,sample=sample,patient=patient,stratum=stratum,gene=gene,n_cells=len(values),mean_log1p=values.mean(),fraction_detected=values.gt(0).mean()))
          if stage=='terminal_DL090':datasets[group,correction]=(umap.to_numpy(),broad,detected.to_numpy())
        z.close()
    qc=pd.DataFrame(qc);expr=pd.DataFrame(expression);comp=pd.DataFrame(composition)
    qc.to_csv(dest/'TCR_stratified_Unknown_QC_by_sample.csv',index=False)
    expr.to_csv(dest/'TCR_stratified_marker_expression_by_sample.csv.gz',index=False)
    comp.to_csv(dest/'lineage_composition_by_sample.csv',index=False)
    def save(fig,name,title,caption,bottom=.09):
        fig.suptitle(title,fontsize=17,fontweight='bold',y=.995)
        fig.text(.012,.012,caption,fontsize=9,va='bottom')
        fig.tight_layout(rect=[0,bottom,1,.965])
        for suffix in ['png','pdf','svg']:fig.savefig(figs/(name+'.'+suffix),bbox_inches='tight')
        figure_rows.append(dict(name=name,title=title,caption=caption,files={s:sha(figs/(name+'.'+s)) for s in ['png','pdf','svg']}))
        plt.close(fig)
    # RNA marker panel, stratified by predicted identity and orthogonal detection.
    strata=[a+' / '+b for a in ['Called T','Unknown','Other called'] for b in ['TCR+','TCR not detected']]
    fig,axs=plt.subplots(2,3,figsize=(19,12))
    for i,group in enumerate(GROUPS):
      for j,correction in enumerate(['NONE','CCA','HARMONY']):
        ax=axs[i,j];q=expr[expr.group.eq(group)&expr.correction.eq(correction)&expr.stage.eq('terminal_DL090')]
        p=q.groupby(['stratum','gene'])[['mean_log1p','fraction_detected']].mean()
        for yi,s in enumerate(strata):
          for xi,g in enumerate(genes):
            if (s,g) not in p.index:continue
            v=p.loc[s,g];scatter=ax.scatter(xi,yi,s=12+150*v.fraction_detected,c=[v.mean_log1p],cmap='viridis',vmin=0,vmax=3.5,edgecolors='none')
        ax.set(xticks=range(len(genes)),xticklabels=genes,yticks=range(len(strata)),yticklabels=strata,title=group+' / '+correction,ylim=(5.6,-.6),xlim=(-.6,len(genes)-.4))
        ax.tick_params(axis='x',labelrotation=70,labelsize=8);ax.grid(axis='y',alpha=.15)
        fig.colorbar(scatter,ax=ax,fraction=.025,pad=.02,label='Mean log1p RNA')
    save(fig,'ptc_RNA_marker_evidence','PTC: RNA marker evidence for called T, Unknown and other cells',
      'Fixed UMAP2 / R HDBSCAN50, AllTissues/mean, terminal DL0.90. Color: mean log1p RNA; dot area increases with detection fraction; means weight available patients equally.\nStrata use productive high-confidence TCR. Missing strata have no dots. RNA modules overlap annotation markers and are descriptive, not independent label truth.')
    # Per-sample/patient composition before and after DL on a fixed path.
    samples=sum(GROUPS.values(),[]);patient_map=ref.groupby('sample').patient.first().to_dict()
    fig,axs=plt.subplots(2,3,figsize=(19,11))
    for i,stage in enumerate(['marker_only_ablation','terminal_DL090']):
      for j,correction in enumerate(['NONE','CCA','HARMONY']):
        ax=axs[i,j];q=comp[comp.correction.eq(correction)&comp.stage.eq(stage)]
        tab=q.pivot(index='sample',columns='broad_lineage',values='fraction').reindex(index=samples,columns=lineage_order).fillna(0)
        base=np.zeros(len(samples))
        for lineage in lineage_order:
            ax.bar(range(8),tab[lineage],bottom=base,color=palette[lineage],label=lineage,width=.8);base+=tab[lineage].to_numpy()
        ax.axvline(3.5,color='black',lw=.8);ax.set(xticks=range(8),xticklabels=[s+'\n'+patient_map[s].replace('Patient','P') for s in samples],ylim=(0,1),ylabel='Cell fraction',title=correction+' / '+stage.replace('_',' '));ax.tick_params(axis='x',labelsize=8)
    handles=[Line2D([],[],marker='s',lw=0,color=palette[l],label=l,markersize=8) for l in lineage_order]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.066),ncol=7,frameon=False,fontsize=9)
    save(fig,'ptc_patient_tissue_composition','PTC: terminal lineage composition by sample, patient and tissue',
      'The same full-RNA AllTissues/mean context and UMAP2 / R HDBSCAN50 partitions are used before and after DL. Vertical lines separate MTN and TUT.\nEach group contains one sample from each patient. Differences are descriptive composition changes, not causal tissue or treatment effects.',bottom=.145)
    # All cells shown: fixed final broad calls and actual TCR detection in each embedding.
    fig,axs=plt.subplots(4,3,figsize=(17,19))
    rng=np.random.default_rng(42)
    for gi,group in enumerate(GROUPS):
      for j,correction in enumerate(['NONE','CCA','HARMONY']):
        u,b,d=datasets[group,correction];order=rng.permutation(len(u));ax=axs[2*gi,j]
        ax.scatter(u[order,0],u[order,1],s=.7,c=[palette[v] for v in b[order]],alpha=.65,rasterized=True,linewidths=0)
        ax.set(title=group+' / '+correction+' / terminal calls',xlabel='UMAP1',ylabel='UMAP2')
        ax=axs[2*gi+1,j];ax.scatter(u[:,0],u[:,1],s=.5,c='#d6d6d6',alpha=.5,rasterized=True,linewidths=0)
        ax.scatter(u[d,0],u[d,1],s=.65,c='#226eae',alpha=.65,rasterized=True,linewidths=0)
        ax.set(title=group+' / '+correction+' / productive TCR+',xlabel='UMAP1',ylabel='UMAP2')
    fig.legend(handles=handles+[Line2D([],[],marker='o',lw=0,color='#226eae',label='TCR+ (detection rows)',markersize=6)],loc='lower center',bbox_to_anchor=(.5,.045),ncol=7,frameon=False,fontsize=9)
    save(fig,'ptc_embedding_TCR_lineage','PTC: pooled embeddings with final calls and TCR detection',
      'All cells are shown. Top row of each group: AllTissues/mean terminal DL0.90 calls from R HDBSCAN50. Lower row: actual productive/high-confidence TCR detection.\nEach embedding has its own coordinates; visual compactness or separation is not a quantitative accuracy measure.',bottom=.10)
    write_json(dest/'manifest.json',dict(status='complete',source_sha256=sha(Path(__file__)),
       label_rules_sha256=sha(Path(__file__).with_name('label_rules.py')),inputs=inputs,figures=figure_rows,
       fixed_context='Pooled UMAP2/R_HDBSCAN50/fullRNA/AllTissues/mean, initial and final0.90',
       marker_genes=genes,n_QC_rows=len(qc),n_expression_rows=len(expr),n_composition_rows=len(comp),
       limitation='No historical doublet score/weights recovered. QC is exact original metadata; RNA evidence descriptive and may overlap marker scoring.',
       outputs={p.name:sha(p) for p in dest.glob('*.csv*')},job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    (dest/'COMPLETE').write_text(sha(dest/'manifest.json')+'\n')
    print('RNA markers, Unknown/TCR QC, patient composition and embeddings complete',flush=True)

if __name__=='__main__':run()
