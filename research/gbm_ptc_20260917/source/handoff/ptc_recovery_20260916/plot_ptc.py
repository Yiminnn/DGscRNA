"""Standalone PTC figures and a complete GBM/PTC decision tree; SLURM only."""
import json
import os
import shutil
from pathlib import Path
from ptc_common import BASE,GROUPS,require_slurm,sha,utc,write_json

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
    from summarize_ptc import FEATURES,ENDPOINTS,effect_values
    preview=os.environ.get('PTC_PLOT_PREVIEW')=='1'
    summary_dir=BASE/('summary_smoke' if preview else 'summary')
    assert (summary_dir/('SMOKE_COMPLETE' if preview else 'COMPLETE')).exists()
    dest=BASE/('figures_preview' if preview else 'figures');dest.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,
        'axes.labelsize':11,'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':180,
        'axes.spines.top':False,'axes.spines.right':False})
    inventory=[]
    def save(fig,name,title,caption):
        fig.suptitle(title,fontsize=17,fontweight='bold',y=.995)
        fig.text(.015,.008,caption,fontsize=9,va='bottom',wrap=True)
        fig.tight_layout(rect=[0,.07,1,.965])
        if preview:fig.text(.98,.97,'PREVIEW: full matrix still pending',ha='right',fontsize=9,color='#a33030')
        for suffix in ['png','pdf','svg']:fig.savefig(dest/(name+'.'+suffix),bbox_inches='tight')
        inventory.append(dict(name=name,title=title,caption=caption,files={s:sha(dest/(name+'.'+s)) for s in ['png','pdf','svg']}))
        plt.close(fig)
    sm=pd.read_csv(summary_dir/'condition_summary.csv.gz',low_memory=False)
    data=pd.read_csv(summary_dir/'all_sample_stage_metrics.csv.gz',low_memory=False)
    t=data[data.stage.eq('terminal_DL090')]
    labels={'cluster_S2_native_ARI':'Clustering ARI vs archived S2',
      'S2_macro_F1_reference_present':'Terminal broad macro-F1 vs S2',
      'productive_strict_recall':'Productive TCR-positive recall',
      'productive_strict_detection_yield':'TCR detection yield in called T cells',
      'productive_strict_apparent_F1':'Apparent TCR binary F1', 'coverage':'Terminal annotation coverage'}
    colors={'matched':'#686868','NONE':'#4C78A8','CCA':'#E68A2E','HARMONY':'#3A9679'}
    # Complete workflow: every experimental intervention is shown at its stage.
    fig,ax=plt.subplots(figsize=(17,13));fig.subplots_adjust(left=.03,right=.97,bottom=.055,top=.95)
    ax.set(xlim=(0,17),ylim=(0,13));ax.axis('off')
    def box(x,y,w,h,text,color='#eef3f7',fs=11):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.07',facecolor=color,edgecolor='#526777',lw=1.1))
        ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=fs)
    def arrow(a,b):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=13,color='#526777',lw=1.3))
    box(.4,11.65,16.2,.8,'Fixed cells and provenance  →  reference-label-free geometry  →  marker scoring + terminal DL/refinement  →  patient-level evaluation','#e0edf4',13)
    box(.5,10.3,7.4,.85,'GBM: 121 samples / 59 patients\nPrimary cohort: 97 samples / 55 patients; exclusions kept separate','#e7f1e9')
    box(9.1,10.3,7.4,.85,'PTC: 92,404 cells / 8 samples / 4 patients\nGate: exact original S2/S3 output replay + full raw-count/QC parity','#f7ede2')
    arrow((4.2,11.57),(4.2,11.23));arrow((12.8,11.57),(12.8,11.23))
    left=[
      ('A1  Geometry features\nAll / VST500 / 1k / 2k / 3k / 5k\nLegacy HVG flavor and marker-union controls',1.0),
      ('A2–A3  Preprocessing and representation\nDirect vs PCA30; seven reductions/no reduction\nUMAP dimensions, legacy bridge, five seeds',1.0),
      ('A4  Clustering\nKMeans / GMM / HDBSCAN15/15\nFixed K23; separate label-free K and density sensitivities',1.0),
      ('A5  Marker scoring\nFull scoring genes held fixed in main contrasts\nSeparate HVG-only scoring and marker-retention ablations',1.0),
      ('A6  Terminal refinement\n15 epochs / threshold 0.90 / fixed all-gene DL input\nMarker-only ablation; no-op / abstention / unavailable explicit',1.0),
      ('GBM output: 39,930 conditions accounted\n35,211 verified terminal outputs\nSingleton sensitivity and fixed-cohort availability bounds',1.0)]
    right=[
      ('P1  Grouping and matched R controls\nMT+N: 48,255 cells | TU+T: 44,149 cells\nEach group has 4 patients; matched single-sample controls',1.0),
      ('P2  Features, representations and correction\nPer sample: six feature levels × seven representations\nUMAP direct vs PCA30; separate PCA30 control\nPooled: NONE | CCA: expression | Harmony: PCA30\nPooled PCA30 and UMAP2 evaluated separately',1.0),
      ('P3  Clustering\nTransfer: KMeans / GMM / Python HDBSCAN15/15\nMatched R: SNN0.5 / HDBSCAN minPts50\nSeparate five-seed Python UMAP/HDBSCAN control',1.0),
      ('P4  Original marker scoring\n17 symbol libraries × none / mean / 0.5 cutoffs\nFull RNA primary; CCA integrated scoring/DL separate',1.0),
      ('P5  Terminal refinement\nOriginal MLP: 10 epochs; fixed RNA2000 input\n0.90 primary / 0.70 sensitivity / marker-only ablation',1.0),
      ('PTC output: 11,696 verified terminal conditions\nSaved-model forward checks; actual/no-op/no-known status\nOne GMM precision rescue; R HDBSCAN index fix verified',1.0)]
    previous=10.3
    for i,((lt,h),(rt,_)) in enumerate(zip(left,right)):
        y=8.9-i*1.25
        box(.5,y,7.4,h,lt,'#f0f6f1',10.5);box(9.1,y,7.4,h,rt,'#fbf3e9',10.5)
        arrow((4.2,previous-.08),(4.2,y+h+.08));arrow((12.8,previous-.08),(12.8,y+h+.08));previous=y
    box(.6,.95,15.8,1.0,'Evaluation after fitting: clustering concordance and terminal annotation reported separately\nGBM: patient-paired effects / availability / held-out HVG choice\nPTC: historical concordance + TCR detection + RNA support + batch/state preservation + held-out panel selection','#e8edf7',11)
    arrow((4.2,2.57),(4.2,2.03));arrow((12.8,2.57),(12.8,2.03))
    fig.suptitle('DG-scRNA: complete workflow and ablation decision tree',fontsize=18,fontweight='bold')
    fig.text(.035,.02,'A1–A6 and P1–P5 mark tested interventions. Reference labels never enter feature selection, embedding, clustering or DL fitting.\nPTC panel choice uses three patients; the fourth is evaluated after selection. Predictions were fitted transductively.',fontsize=10)
    if preview:fig.text(.98,.965,'PREVIEW: full matrix still pending',ha='right',fontsize=9,color='#a33030')
    for suffix in ['png','pdf','svg']:fig.savefig(dest/('decision_tree_complete.'+suffix),bbox_inches='tight')
    inventory.append(dict(name='decision_tree_complete',title='Complete workflow and ablation decision tree',
        caption='All GBM and PTC stages, tested interventions, terminal-output rules and validation layers.',
        files={s:sha(dest/('decision_tree_complete.'+s)) for s in ['png','pdf','svg']}));plt.close(fig)
    # Original counts, QC and TCR definitions.
    qc=pd.read_csv(summary_dir/'raw_counts_QC_parity.csv')
    fig,ax=plt.subplots(1,2,figsize=(15,6))
    x=np.arange(len(qc));w=.25
    for i,(col,lab,c) in enumerate([('n_raw_cells','Raw 10x','#aab6c5'),('n_after_raw_reference_QC','QC on full raw genes','#6093bd'),('n_archived_S2','Archived retained','#277865')]):
        ax[0].bar(x+(i-1)*w,qc[col],w,label=lab,color=c)
    ax[0].set(xticks=x,xticklabels=qc['sample'],ylabel='Cells',ylim=(0,qc.n_raw_cells.max()*1.3),title='Exact retained counts and QC; original fixed cells');ax[0].legend(fontsize=9)
    tc=pd.read_csv(BASE/'archived_comparators/historical_TCR_endpoint_replay_corrected_names.csv')
    rr=tc[tc.method.eq('S2_terminal')&tc['mapping'].eq('strict_T')&tc['sample'].eq('ALL')]
    definitions=['TCR_any_filtered_contig','TCR_cell_high_confidence_productive_TCR','TCR_paired_productive_TRA_TRB','TCR_S3_supplied']
    vals=rr.set_index('TCR_definition').loc[definitions,'n_detected_TCR']
    bars=ax[1].bar(['Any contig\nS2 definition','Productive / HC\nprimary','Paired TRA/TRB','S3 supplied\nversion sensitivity'],vals,color=['#6093bd','#277865','#75ae9d','#cf9270'])
    ax[1].bar_label(bars,fmt='%d',padding=4);ax[1].set(ylabel='TCR-detected cells',ylim=(0,42000),title='TCR detection definitions remain separate')
    save(fig,'ptc_source_QC_audit','PTC recovery: source parity before new experiments',
      'All 33,694 raw genes match in all 92,404 retained cells after 15 Seurat symbol substitutions. Archived outputs reproduce exactly; historical model weights are absent.\nQC uses full raw-gene metadata; a gene-filter-first candidate loses 3 additional cells. Count compatibility does not recover the original stochastic doublet model.')
    # Complete transfer factorial.
    f=pd.read_csv(summary_dir/'primary_factorial_patient_means.csv');f=f[f.analysis_scope.eq('ALL_8_SAMPLES')]
    f['method']=f.dr+' / '+f.clusterer
    order=[d+' / '+c for d in ['PCA','FA','ICA','ISOMAP','UMAP','TSNE','none'] for c in ['KMeans','GMM','HDBSCAN']]
    # The implementation records Isomap with mixed case; retain exact names.
    order=[d+' / '+c for d in ['PCA','FA','ICA',next(z for z in f.dr.unique() if z.lower()=='isomap'),'UMAP','TSNE','none'] for c in ['KMeans','GMM','HDBSCAN']]
    fig,axs=plt.subplots(2,2,figsize=(19,12))
    for ax,metric in zip(axs.flat,[ENDPOINTS[0],ENDPOINTS[1],ENDPOINTS[2],ENDPOINTS[5]]):
        table=f.pivot(index='feature',columns='method',values=metric).reindex(index=FEATURES,columns=order)
        im=ax.imshow(table,aspect='auto',vmin=-.05 if metric==ENDPOINTS[0] else 0,vmax=1,cmap='viridis')
        ax.set(xticks=np.arange(21),xticklabels=order,yticks=np.arange(6),yticklabels=['All','500','1,000','2,000','3,000','5,000'],title=labels[metric],ylabel='Geometry features')
        ax.tick_params(axis='x',labelrotation=65,labelsize=8)
        fig.colorbar(im,ax=ax,fraction=.025,pad=.02)
    save(fig,'ptc_complete_factorial','PTC: complete HVG × representation × clustering comparison',
      'Terminal DL/refinement 0.90; fixed AllTissues/mean, full RNA scoring and group-selected RNA2000 classifier input. Patient-balanced means across 8 samples / 4 patients.\nHistorical concordance and TCR-positive recall are distinct endpoints. The N-2/HVG2000/t-SNE/GMM precision recovery is flagged in the tables.')
    # HVG dose/path response with exhaustive patient-bootstrap intervals.
    primary=t[t['mode'].eq('single')&t.library.eq('CellMarker_AllTissues')&t.cutoff.eq('mean')&t.seed.eq(42)&t.clusterer.eq('HDBSCAN')]
    fig,axs=plt.subplots(2,2,figsize=(14,10))
    for ax,metric in zip(axs.flat,[ENDPOINTS[0],ENDPOINTS[1],ENDPOINTS[2],ENDPOINTS[5]]):
        for dr,dim,inp,name,c,offset in [('UMAP',2,'direct','Direct UMAP2','#4C78A8',-.08),('UMAP',2,'pca30','PCA30 → UMAP2','#E68A2E',0),('PCA',30,'direct','PCA30','#3A9679',.08)]:
            q=primary[primary.dr.eq(dr)&primary.dim.eq(dim)&primary.input_space.eq(inp)]
            xx=[];yy=[];err=[]
            for j,feature in enumerate(FEATURES):
                v=q[q.feature.eq(feature)].groupby('patient')[metric].mean()
                if len(v):
                    e=effect_values(v);xx.append(j+offset);yy.append(e['delta_mean']);err.append([e['delta_mean']-e['ci_lower'],e['ci_upper']-e['delta_mean']])
            ax.errorbar(xx,yy,yerr=np.asarray(err).T,fmt='o-',capsize=3,color=c,label=name,alpha=.9)
        ax.set(xticks=range(6),xticklabels=['All','500','1k','2k','3k','5k'],xlabel='Geometry features',ylabel=labels[metric]);ax.grid(alpha=.15);ax.legend(fontsize=9)
    save(fig,'ptc_HVG_PCA_response','PTC: HVG effects depend on the preprocessing path',
      'Points are four-patient means; intervals exhaust all 256 size-four patient-bootstrap resamples and are descriptive. PCA paths were tested at all genes and HVG2000 only.\nThe paired sign-flip test has only 16 permutations: minimum two-sided p = 0.125. Neither cells nor seeds increase the biological sample size.')
    # Matched R single/pooling/correction comparisons, fixed representative path.
    q=t[t['mode'].isin(['matched','pooled'])&t.scoring_assay.eq('RNA')&t.library.eq('CellMarker_AllTissues')&t.cutoff.eq('mean')&t.space.eq('UMAP2')&t.clusterer.eq('HDBSCAN_R')].copy()
    q['workflow']=np.where(q['mode'].eq('matched'),'matched',q.correction)
    fig,axs=plt.subplots(2,4,figsize=(18,9))
    workflow=['matched','NONE','CCA','HARMONY']
    for row,group in enumerate(GROUPS):
      for ax,metric in zip(axs[row],[ENDPOINTS[1],ENDPOINTS[4],ENDPOINTS[2],ENDPOINTS[5]]):
        z=q[q.group.eq(group)]
        for patient,p in z.groupby('patient'):
            vals=p.set_index('workflow').loc[workflow,metric];ax.plot(range(4),vals,'o-',color='#929aa3',alpha=.65,lw=1,ms=4)
        means=z.groupby('workflow')[metric].mean().reindex(workflow)
        ax.plot(range(4),means,'o-',color='#1d5673',lw=2.5,ms=6)
        ax.set(xticks=range(4),xticklabels=['Single','Pooled\nNONE','CCA','Harmony'],ylim=(-.02,1.03),title=labels[metric],ylabel=group);ax.grid(axis='y',alpha=.15)
    save(fig,'ptc_matched_pooling_batch','PTC: separate pooling from correction using matched R controls',
      'Fixed PCA30 → UMAP2 / R HDBSCAN minPts50, AllTissues/mean, full RNA scoring and constant RNA2000 DL input. Gray lines are the four paired patients per group.\nAll PCA30/SNN combinations are also reported in the tables. Sample, tissue and patient are confounded; greater batch mixing alone does not establish correct biology.')
    # Full 17-library x 3-cutoff roster for a fixed, preselected R path.
    for metric,name in [(ENDPOINTS[4],'ptc_marker_roster_TCR'),(ENDPOINTS[5],'ptc_marker_roster_coverage')]:
        fig,axs=plt.subplots(2,3,figsize=(15,15))
        libraries=sorted(t.library.unique())
        for i,group in enumerate(GROUPS):
          for j,correction in enumerate(['NONE','CCA','HARMONY']):
            ax=axs[i,j];z=t[t['mode'].eq('pooled')&t.group.eq(group)&t.correction.eq(correction)&t.scoring_assay.eq('RNA')&t.space.eq('UMAP2')&t.clusterer.eq('HDBSCAN_R')]
            values=z.groupby(['library','cutoff'])[metric].mean().unstack().reindex(index=libraries,columns=['none','mean','0.5'])
            im=ax.imshow(values,aspect='auto',vmin=0,vmax=1,cmap='viridis')
            ax.set(xticks=range(3),xticklabels=['None','Mean','0.5'],yticks=range(len(libraries)),yticklabels=[s.replace('CellMarker_','CM: ').replace('Pubmed_34663816','GSE184362-derived') for s in libraries],title=group+' / '+correction,xlabel='Density cutoff')
            for y in range(len(libraries)):
                for x in range(3):ax.text(x,y,f'{values.iloc[y,x]:.2f}',ha='center',va='center',fontsize=9,color='white' if values.iloc[y,x]<.65 else 'black')
            fig.colorbar(im,ax=ax,fraction=.03,pad=.02)
        save(fig,name,'PTC original 17-library roster: '+labels[metric],
          'Each cell is the mean across four patients, terminal DL/refinement 0.90. Fixed UMAP2 / R HDBSCAN minPts50, full RNA scoring, RNA2000 DL input.\nThis descriptive grid is not a selection-adjusted winner ranking. All libraries, cutoffs, PCA30 and SNN branches are retained in the complete condition tables.')
    # Refinement and threshold plots: condition/sample matched, no aggregate accuracy claim.
    fig,axs=plt.subplots(2,3,figsize=(15,10))
    indices=['condition_id','sample']
    p0=data[data.stage.eq('marker_only_ablation')].set_index(indices)
    p9=t.set_index(indices);p7=data[data.stage.eq('terminal_DL070_sensitivity')].set_index(indices)
    for col,metric in enumerate([ENDPOINTS[1],ENDPOINTS[4],ENDPOINTS[5]]):
      for row,(a,b,xlab,ylab) in enumerate([(p0,p9,'Marker-only ablation','Terminal DL 0.90'),(p9,p7,'Terminal DL 0.90','Terminal DL 0.70')]):
        ax=axs[row,col];a=a.reindex(b.index)
        im=ax.hexbin(a[metric],b[metric],gridsize=40,extent=(0,1,0,1),mincnt=1,bins='log',cmap='Blues')
        ax.plot([0,1],[0,1],'--',color='#8b8b8b',lw=1);ax.set(xlabel=xlab,ylabel=ylab,title=labels[metric],xlim=(0,1),ylim=(0,1))
        fig.colorbar(im,ax=ax,label='Sample-condition count (log scale)',fraction=.045)
    save(fig,'ptc_refinement_threshold','PTC: terminal refinement and confidence-threshold ablations',
      'Identical conditions and cells are paired. Marker-only outputs are ablations; every headline DG-scRNA result is terminal. No-op and no-known-label states remain explicit.\nChanging the threshold reuses the same saved probabilities. Plot counts are sample-conditions, not independent biological replicates.')
    # Balanced batch mixing with local and transcript-state preservation.
    bio=pd.read_csv(BASE/'batch_biology/shared_lineage_metrics.csv');bio=bio[bio.scope.eq('productive_TCR_positive')]
    states=['Proliferation','Interferon','Stress','CD4','CD8','Treg']
    fig,axs=plt.subplots(2,2,figsize=(15,10))
    for row,group in enumerate(GROUPS):
        q=bio[bio.group.eq(group)];ax=axs[row,0]
        for r in q.itertuples():
            ax.scatter(r.mean_batch_entropy,r.within_sample_neighbor_retention_vs_NONE,c=colors[r.correction],marker='o' if r.space=='PCA30' else '^',s=90)
            left=r.correction=='NONE' and r.space=='PCA30'
            ax.annotate(r.correction+' / '+r.space,(r.mean_batch_entropy,r.within_sample_neighbor_retention_vs_NONE),
                        xytext=(-6,7) if left else (6,7),ha='right' if left else 'left',textcoords='offset points',fontsize=8)
        ax.set(xlabel='Balanced local batch entropy (higher = more mixing)',ylabel='Within-sample neighbor retention vs NONE',title=group+' / productive TCR-positive cells',ylim=(0,1.1),xlim=(0,1.12));ax.grid(alpha=.2)
        table=q.set_index(['correction','space'])[[s+'_mean_log1p_neighbor_spearman' for s in states]]
        im=axs[row,1].imshow(table,aspect='auto',vmin=-1,vmax=1,cmap='coolwarm')
        axs[row,1].set(xticks=range(6),xticklabels=states,yticks=range(len(table)),yticklabels=[' / '.join(i) for i in table.index],title='RNA-state vs neighbor-average Spearman correlation')
        fig.colorbar(im,ax=axs[row,1],fraction=.04)
    save(fig,'ptc_batch_biology','PTC: batch mixing accompanied by biology-preservation diagnostics',
      'Fixed balanced cells: up to 500 per sample in each shared lineage; 30-neighbor mixing and 15-neighbor within-sample retention. All eligible shared S2 lineages are in the tables.\nRNA modules are descriptive and may overlap annotation markers. Neither these correlations nor TCR absence supplies independent ground-truth cell identity.')
    # Held-out selection across all predeclared pooled paths is available in tables.
    hs=pd.read_csv(summary_dir/'patient_heldout_panel_selection.csv')
    hs=hs[hs.space.eq('UMAP2')&hs.clusterer.eq('HDBSCAN_R')]
    policies=['heldout_selected','CellMarker_AllTissues_mean','CellMarker_Thyroid_mean','Pubmed_34663816_mean','HPA_allThyroid_mean']
    fig,axs=plt.subplots(2,3,figsize=(16,9))
    for row,group in enumerate(GROUPS):
      for col,correction in enumerate(['NONE','CCA','HARMONY']):
        ax=axs[row,col];q=hs[hs.group.eq(group)&hs.correction.eq(correction)]
        for j,metric in enumerate([ENDPOINTS[4],ENDPOINTS[2],ENDPOINTS[5]]):
            means=q.groupby('policy')[metric].mean().reindex(policies)
            ax.plot(range(5),means,'o-',label=labels[metric],lw=1.7,ms=5)
        ax.set(xticks=range(5),xticklabels=['Held-out\nselection','AllTissues','Thyroid','GSE184362\nmarkers','HPA thyroid'],title=group+' / '+correction,ylim=(-.02,1.03));ax.tick_params(axis='x',labelsize=8)
        if row==0 and col==0:ax.legend(fontsize=8)
    save(fig,'ptc_heldout_context','PTC: held-out patient selection of marker library and cutoff',
      'Other three patients choose the panel/cutoff by apparent productive-TCR binary F1; the held-out patient alone supplies each evaluation value. Fixed comparator panels use mean cutoff.\nSelection operates on transductively fitted predictions; expression embeddings and DL models were not refitted without the held-out patient. This is not fully inductive validation.')
    # Archive competitors retain their old workflow and source definitions.
    methods=['S2_initial','S2_terminal','S3_terminal','SCINA_archived','scCATCH_archived','scType_archived','SignacX_archived']
    fig,axs=plt.subplots(1,2,figsize=(16,7))
    for ax,metric,title in zip(axs,['TCR_positive_recall','apparent_TCR_binary_F1'],['Productive TCR-positive recall','Apparent productive-TCR binary F1']):
      for mapping,c,offset in [('strict_T','#4C78A8',-.14),('legacy_broad_sensitivity','#E68A2E',.14)]:
        q=tc[tc['sample'].eq('ALL')&tc['mapping'].eq(mapping)&tc.TCR_definition.eq('TCR_cell_high_confidence_productive_TCR')].set_index('method')
        ax.bar(np.arange(len(methods))+offset,q.loc[methods,metric],.28,label=mapping.replace('_',' '),color=c)
      ax.set(xticks=range(len(methods)),xticklabels=[m.replace('_archived','').replace('_','\n') for m in methods],title=title,ylim=(0,1.05));ax.legend(fontsize=9)
    save(fig,'ptc_archived_comparator_audit','PTC archived methods: original conditions and name-mapping sensitivity',
      'Historical outputs only, pooled over 92,404 cells; this is not a matched ranking against the new two-group reconstruction. Strict T excludes NK/NKT and ambiguous lymphoid names.\nCompact names and all 623 SCINA encoding differences are normalized through the frozen ontology. Literal source binary flags are separate; unresolved S3 TCR/DG-flag discrepancies remain documented.')
    assert (BASE/'biology_detail/COMPLETE').exists()
    detail=json.loads((BASE/'biology_detail/manifest.json').read_text())
    if preview:
        for item in detail['figures']:
            for suffix in item['files']:shutil.copy2(BASE/'figures'/(item['name']+'.'+suffix),dest/(item['name']+'.'+suffix))
    inventory.extend(detail['figures'])
    write_json(dest/'figure_manifest.json',dict(status='preliminary_partial_diagnostic' if preview else 'complete',figures=inventory,n_figures=len(inventory),
        source_sha256=sha(Path(__file__)),summary_manifest_sha256=sha(summary_dir/'manifest.json'),
        completed_at=utc(),job=os.environ['SLURM_JOB_ID']))
    (dest/('SMOKE_COMPLETE' if preview else 'COMPLETE')).write_text(sha(dest/'figure_manifest.json')+'\n')
    print('Rendered',len(inventory),'figures in PNG/PDF/SVG',flush=True)

if __name__=='__main__':run()
