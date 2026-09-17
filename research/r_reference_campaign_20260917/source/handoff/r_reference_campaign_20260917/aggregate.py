"""Aggregate the complete grid without concealing failed/no-op or abstaining arms."""
import json,os,sys,hashlib,time
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
    dest=OUT/'summary';dest.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
    expected=[(u,'HCL' if u.startswith('HCL__') else ('colorectal' if u=='colorectal_CCA28_ge31' else u),OUT/'benchmark'/u/'reference_CCA2000') for u in (OUT/'benchmark_units.txt').read_text().splitlines()]
    specs=json.loads((ROOT/'handoff/r_reference_campaign_20260917/ptc_ablation_conditions.json').read_text())
    expected += [(s['unit'],'PTC',OUT/'PTC_ablation'/s['unit']) for s in specs]
    expected += [(f'PTC_{g}_GEOMETRY{h}_FIXED_CCAall_DL2000','PTC',OUT/'PTC_ablation'/f'PTC_{g}_GEOMETRY{h}_FIXED_CCAall_DL2000') for g in ['NMT','TTU'] for h in ['500','1000','2000','3000','5000','all']]
    expected += [('PTC_archived_CCA2000','PTC',OUT/'PTC_archived_CCA2000')]
    inventories=[];metrics=[];clusters=[];statuses=[];donors=[];abstentions=[]
    grid_path=OUT/'verification/frozen_grid_audit.json'
    grid=json.loads(grid_path.read_text()) if grid_path.exists() else {}
    grid_records={r['unit']:r for r in grid.get('records',[])}
    for unit,dataset,prep in expected:
        pm=json.loads((prep/'prepare_manifest.json').read_text()) if (prep/'prepare_manifest.json').exists() else {}
        if unit=='PTC_archived_CCA2000':pm={'n_cells':92404,'correction':'archived_all8_CCA2000','features':{'anchor':2000,'geometry':2000,'scoring':2000,'DL':2000}}
        manifests=list(prep.glob('*/score_manifest.json'))
        n_expected=sum(len(json.loads(p.read_text())['arms']) for p in manifests)
        n_terminal=len(list(prep.glob('*/terminal/*/TERMINAL_COMPLETE')))
        complete=(prep/'evaluation/COMPLETE').exists()
        evaluation_manifest=prep/'evaluation/evaluation_manifest.json'
        grid_fresh=(unit in grid_records and evaluation_manifest.exists() and
            grid_records[unit]['evaluation_manifest_sha256']==hashlib.sha256(evaluation_manifest.read_bytes()).hexdigest())
        inv=dict(dataset=dataset,unit=unit,directory=str(prep),n_cells=pm.get('n_cells'),
            correction=pm.get('correction'),prepared=bool(pm),scored_routes=len(manifests),
            expected_terminal_from_scored_routes=n_expected,terminal_complete=n_terminal,evaluated=complete,
            plotted=(prep/'figures/FIGURES_COMPLETE').exists(),audited=(prep/'verification/AUDIT_COMPLETE').exists(),
            frozen_grid_audited=grid_fresh,
            abstention_audited=dataset!='PTC' or (prep/'evaluation/ABSTENTION_COMPLETE').exists())
        inv.update({f'features_{k}':v for k,v in pm.get('features',{}).items()})
        inventories.append(inv)
        if not complete:continue
        metrics.append(pd.read_csv(prep/'evaluation/metrics.csv.gz'))
        if (prep/'evaluation/ABSTENTION_COMPLETE').exists():abstentions.append(pd.read_csv(prep/'evaluation/abstention_aware_metrics.csv.gz'))
        for name,collection in [('clustering_metrics.csv',clusters),('terminal_statuses.csv',statuses),('per_donor.csv.gz',donors)]:
            if (prep/'evaluation'/name).exists():
                frame=pd.read_csv(prep/'evaluation'/name);frame['unit']=unit;frame['dataset']=dataset;collection.append(frame)
    inv=pd.DataFrame(inventories);inv.to_csv(dest/'completion_inventory.csv',index=False)
    allm=pd.concat(metrics,ignore_index=True) if metrics else pd.DataFrame()
    if abstentions:
        ab=pd.concat(abstentions,ignore_index=True);ab.to_csv(dest/'PTC_abstention_aware_metrics.csv.gz',index=False)
        keys=['unit','route','arm_id','stage','scope','endpoint']
        extra=[c for c in ab.columns if c not in allm.columns]
        allm=allm.merge(ab[keys+extra],on=keys,how='left',validate='one_to_one')
    allm.to_csv(dest/'all_annotation_metrics.csv.gz',index=False)
    cl=pd.concat(clusters,ignore_index=True) if clusters else pd.DataFrame();cl.to_csv(dest/'all_clustering_metrics.csv',index=False)
    st=pd.concat(statuses,ignore_index=True) if statuses else pd.DataFrame();st.to_csv(dest/'all_terminal_statuses.csv',index=False)
    if not len(allm):return
    bench=allm[allm.dataset.ne('PTC')]
    fixed=bench[bench.library.eq('CM2_primary_normal') & bench.cutoff.eq('mean')]
    fixed.to_csv(dest/'benchmark_fixed_primary_normal_mean.csv',index=False)
    final=bench[bench.stage.eq('final090')]
    maximum=[]
    for (dataset,unit,endpoint),group in final.groupby(['dataset','unit','endpoint']):
        best=group.sort_values(['macro_F1','accuracy','unknown_fraction','route','library','cutoff'],ascending=[False,False,True,True,True,True]).iloc[0]
        maximum.append(dict(dataset=dataset,unit=unit,endpoint=endpoint,n_tested=len(group),
          median_macro_F1=group.macro_F1.median(),descriptive_maximum_macro_F1=best.macro_F1,
          maximum_route=best.route,maximum_library=best.library,maximum_cutoff=best.cutoff,
          accuracy_at_maximum=best.accuracy,unknown_at_maximum=best.unknown_fraction,
          interpretation='Descriptive best on these same evaluation labels; not an independently selected operating point'))
    pd.DataFrame(maximum).to_csv(dest/'benchmark_descriptive_maxima.csv',index=False)
    # HCL summaries explicitly average tissue-specific metrics, never masquerading as a pooled fit.
    hcl=[]
    for keys,group in bench[bench.dataset.eq('HCL')].groupby(['library','route','cutoff','stage','endpoint']):
        row=dict(zip(['library','route','cutoff','stage','endpoint'],keys))
        row.update(n_tissue_units=group.unit.nunique(),n_cells=group.n_cells.sum(),
          equal_tissue_mean_macro_F1=group.macro_F1.mean(),
          cell_weighted_mean_tissue_macro_F1=np.average(group.macro_F1,weights=group.n_cells),
          pooled_accuracy=np.average(group.accuracy,weights=group.n_cells),
          pooled_unknown_fraction=np.average(group.unknown_fraction,weights=group.n_cells))
        hcl.append(row)
    pd.DataFrame(hcl).to_csv(dest/'HCL_tissue_stratified_summary.csv',index=False)
    # A donor-held-out label-selection sensitivity. All cells enter unsupervised fitting;
    # held-out donor labels alone are excluded from choosing the marker/route/cutoff.
    selections=[]
    if donors:
        d=pd.concat(donors,ignore_index=True);d.to_csv(dest/'all_per_donor.csv.gz',index=False)
        d=d[d.stage.eq('final090')]
        keys=['route','library','cutoff']
        for unit,group in d.groupby('unit'):
            if group.donor.nunique()<3:continue
            for donor in sorted(group.donor.unique()):
                train=group[group.donor.ne(donor)]
                ranks=train.groupby(keys,as_index=False).macro_F1.mean().sort_values(['macro_F1',*keys],ascending=[False,True,True,True])
                best=ranks.iloc[0]
                held=group[group.donor.eq(donor)]
                for key in keys:held=held[held[key].eq(best[key])]
                assert len(held)==1
                record=held.iloc[0].to_dict();record['selection_training_donors_mean_macro_F1']=best.macro_F1
                record['selection_scope']='Held-out annotation labels; transductive fixed cohort fit, no donor-independent refit'
                selections.append(record)
        pd.DataFrame(selections).to_csv(dest/'donor_heldout_label_selection.csv',index=False)
    ptc=allm[allm.dataset.eq('PTC')].copy()
    primary=ptc[((ptc.scope.eq('NMT') & ptc.library.eq('CellMarker_Thyroid') & ptc.cutoff.eq('none')) |
                 (ptc.scope.eq('TTU') & ptc.library.eq('Pubmed_34663816') & ptc.cutoff.eq('mean')))]
    primary.to_csv(dest/'PTC_fixed_historical_marker_contexts.csv',index=False)
    primary[primary.stage.eq('final090') & primary.endpoint.eq('strict_T_name_rule')].to_csv(dest/'PTC_terminal_primary_contexts_strict_T.csv',index=False)
    ptc_max=[]
    pq=ptc[ptc.stage.eq('final090') & ptc.endpoint.eq('strict_T_name_rule') & ptc.scope.isin(['NMT','TTU'])].copy()
    def family(unit):
        if unit=='PTC_archived_CCA2000':return 'archived_geometry_refit'
        if unit=='PTC_ALL8_CCA2000':return 'fresh_all8_refit'
        if '_GEOMETRY' in unit:return 'isolated_geometry'
        if '_RNA2000' in unit:return 'RNA_NONE_HARMONY'
        return 'joint_CCA_budget'
    pq['family']=pq.unit.map(family)
    for (group,mode),frame in pq.groupby(['scope','family']):
        for metric in ['F1_T','macro_F1_T_nonT_unknown_as_error','accuracy_unknown_as_error']:
            if metric not in frame or frame[metric].isna().all():continue
            winner=frame.sort_values([metric,'unknown_fraction','unit','route','library','cutoff'],ascending=[False,True,True,True,True,True]).iloc[0]
            ptc_max.append(dict(group=group,family=mode,metric=metric,value=winner[metric],unit=winner.unit,
                route=winner.route,library=winner.library,cutoff=winner.cutoff,unknown_fraction=winner.unknown_fraction,
                evaluated_conditions=len(frame),interpretation='Descriptive same-label maximum within this comparison family; not a universal optimum'))
    pd.DataFrame(ptc_max).to_csv(dest/'PTC_descriptive_maxima_by_family.csv',index=False)
    numeric=['macro_F1','weighted_F1','F1_T','F1_nonT','AUC_binary','accuracy','unknown_fraction','saved_native_concordance','saved_broad_lineage_concordance',
             'accuracy_unknown_as_error','F1_nonT_unknown_as_error','macro_F1_T_nonT_unknown_as_error','coverage']
    keys=['dataset','unit','route','library','cutoff','endpoint','scope']
    before=allm[allm.stage.eq('initial')];after=allm[allm.stage.eq('final090')]
    changes=after.merge(before,on=keys,suffixes=('_final','_initial'),validate='one_to_one')
    for c in numeric:
        if c+'_final' in changes:changes['delta_'+c]=changes[c+'_final']-changes[c+'_initial']
    changes.to_csv(dest/'DL_vs_initial_paired_changes.csv.gz',index=False)
    ptc['route_canonical']=ptc.route.replace({'seurat_clusters':'PCA30_SNN','seurat.UMAP_clusters':'UMAP2_SNN','hdbscan_clusters':'PCA30_HDBSCAN_R','hdbscan.UMAP_clusters':'UMAP2_HDBSCAN_R'})
    paired=[]
    pairkeys=['scope','route_canonical','library','cutoff','stage','endpoint']
    for group in ['NMT','TTU']:
        scope=ptc[ptc.scope.eq(group)]
        for family,baseline,conditions in [
          ('joint_CCA_budget',f'PTC_{group}_CCA2000',[f'PTC_{group}_CCA{h}' for h in ['500','1000','2000','3000','5000','all']]),
          ('geometry_only',f'PTC_{group}_GEOMETRY2000_FIXED_CCAall_DL2000',[f'PTC_{group}_GEOMETRY{h}_FIXED_CCAall_DL2000' for h in ['500','1000','2000','3000','5000','all']]),
          ('matched_RNA_geometry_correction',f'PTC_{group}_NONE_RNA2000',[f'PTC_{group}_HARMONY_RNA2000']),
          ('combined_CCA_vs_RNA',f'PTC_{group}_CCA2000',[f'PTC_{group}_NONE_RNA2000',f'PTC_{group}_HARMONY_RNA2000']),
          ('integration_scope_and_DL_pool','PTC_ALL8_CCA2000',[f'PTC_{group}_CCA2000']),
          ('fresh_vs_archived_all8','PTC_archived_CCA2000',['PTC_ALL8_CCA2000'])]:
            ref=scope[scope.unit.eq(baseline)]
            if ref.empty:continue
            for condition in conditions:
                q=scope[scope.unit.eq(condition)]
                if q.empty:continue
                p=q.merge(ref,on=pairkeys,suffixes=('_condition','_baseline'),validate='one_to_one')
                assert (p.n_cells_condition==p.n_cells_baseline).all()
                p['contrast_family']=family
                for c in numeric:
                    if c+'_condition' in p:p['delta_'+c]=p[c+'_condition']-p[c+'_baseline']
                paired.append(p)
    if paired:pd.concat(paired,ignore_index=True).to_csv(dest/'PTC_paired_ablation_changes.csv.gz',index=False)
    palette={'PCA30_SNN':'#0072B2','PCA30_HDBSCAN_R':'#D55E00','UMAP2_SNN':'#009E73','UMAP2_HDBSCAN_R':'#CC79A7'}
    symbols={'PCA30_SNN':'o','PCA30_HDBSCAN_R':'s','UMAP2_SNN':'^','UMAP2_HDBSCAN_R':'D'}
    short={'PCA30_SNN':'PCA + SNN','PCA30_HDBSCAN_R':'PCA + HDBSCAN','UMAP2_SNN':'UMAP + SNN','UMAP2_HDBSCAN_R':'UMAP + HDBSCAN'}
    def save(fig,name):
        for ext in ['png','pdf']:fig.savefig(dest/(name+'.'+ext),dpi=300,bbox_inches='tight',facecolor='white')
        plt.close(fig)
    fp=fixed[fixed.stage.eq('final090') & fixed.endpoint.eq('common_lineage') & fixed.dataset.ne('HCL')]
    if len(fp):
        piv=fp.pivot(index='unit',columns='route',values='macro_F1').reindex(columns=list(palette))
        hf=pd.DataFrame(hcl)
        if len(hf):
            hs=hf[hf.library.eq('CM2_primary_normal') & hf.cutoff.eq('mean') & hf.stage.eq('final090') & hf.endpoint.eq('common_lineage')]
            if len(hs)==4:piv.loc[f'HCL (mean of {int(hs.n_tissue_units.min())} tissue groups)']=hs.set_index('route').loc[piv.columns,'equal_tissue_mean_macro_F1'].to_numpy()
        fig,ax=plt.subplots(figsize=(7.2,max(3.5,.32*len(piv)+1.5)))
        im=ax.imshow(piv.to_numpy(),vmin=0,vmax=1,cmap='cividis',aspect='auto')
        ax.set_xticks(range(4));ax.set_xticklabels([short[x] for x in piv.columns]);ax.set_yticks(range(len(piv)));ax.set_yticklabels(piv.index)
        for i in range(len(piv)):
            for j in range(4):
                v=piv.iloc[i,j]
                if np.isfinite(v):ax.text(j,i,f'{v:.3f}',ha='center',va='center',color='white' if v<.5 else 'black')
        fig.colorbar(im,ax=ax,label='Common-lineage macro F1',fraction=.04,pad=.025)
        ax.set_title('Reference workflow | fixed primary-normal CellMarker / mean\nTerminal DL; no best-on-test marker selection',fontsize=10)
        fig.tight_layout();save(fig,'reviewer_fixed_context_comparison')
    for mode in ['joint_CCA_budget','geometry_only']:
        fig,axs=plt.subplots(3,2,figsize=(7.2,7.7),sharex=True)
        any_data=False
        for col,group in enumerate(['NMT','TTU']):
            for route,color in palette.items():
                rows=[]
                for budget in ['500','1000','2000','3000','5000','all']:
                    unit=f'PTC_{group}_CCA{budget}' if mode=='joint_CCA_budget' else f'PTC_{group}_GEOMETRY{budget}_FIXED_CCAall_DL2000'
                    q=primary[primary.unit.eq(unit)&primary.route.eq(route)&primary.stage.eq('final090')&primary.endpoint.eq('strict_T_name_rule')]
                    if len(q):rows.append((budget,q.iloc[0]));any_data=True
                if not rows:continue
                xx=[['500','1000','2000','3000','5000','all'].index(b) for b,_ in rows]
                for row,metric in enumerate(['F1_T','saved_broad_lineage_concordance','unknown_fraction']):
                    axs[row,col].plot(xx,[r[metric] for _,r in rows],marker=symbols[route],ms=3,lw=1.2,color=color,label=short[route])
            axs[0,col].set_title(group+(' | Thyroid / none' if group=='NMT' else ' | Pubmed / mean'),fontsize=9)
            for row in range(3):
                axs[row,col].set_ylim(0,1);axs[row,col].grid(axis='y',alpha=.18);axs[row,col].spines[['top','right']].set_visible(False)
                axs[row,col].set_xticks(range(6));axs[row,col].set_xticklabels(['500','1k','2k','3k','5k','All'])
            axs[2,col].set_xlabel('Requested geometry gene budget' if mode=='geometry_only' else 'Requested CCA feature budget')
        for row,label in enumerate(['TCR agreement: strict T F1','Archived lineage concordance','Unknown fraction']):axs[row,0].set_ylabel(label)
        handles,labels=axs[0,0].get_legend_handles_labels()
        if not handles:handles,labels=axs[0,1].get_legend_handles_labels()
        fig.legend(handles,labels,loc='lower center',ncol=4,fontsize=7,bbox_to_anchor=(.5,-.012))
        title='PTC | joint anchor / geometry / score / DL feature-budget intervention' if mode=='joint_CCA_budget' else 'PTC | geometry-only intervention on fixed all-gene CCA expression'
        fig.suptitle(title,fontsize=10);fig.tight_layout(rect=(0,.025,1,.96))
        if any_data:save(fig,'PTC_'+mode)
        else:plt.close(fig)
    complete=bool(inv.evaluated.all() and inv.plotted.all() and inv.audited.all() and inv.abstention_audited.all()
        and inv.frozen_grid_audited.all() and grid.get('status')=='passed' and grid.get('units')==len(inv))
    report=dict(status='complete' if complete else 'in_progress',timestamp_epoch=time.time(),job=os.environ['SLURM_JOB_ID'],
       expected_analysis_units=len(inv),evaluated_units=int(inv.evaluated.sum()),plotted_units=int(inv.plotted.sum()),audited_units=int(inv.audited.sum()),
       terminal_annotation_conditions_evaluated=len(st),clustering_conditions_evaluated=len(cl),
       expected_reviewer_datasets=11,expected_HCL_tissue_units=59,expected_PTC_units=30,
       missing_evaluation=inv.loc[~inv.evaluated,'unit'].tolist(),missing_figures=inv.loc[~inv.plotted,'unit'].tolist(),
       terminal_status_counts=st.dl_status.value_counts().to_dict() if len(st) else {},
       no_global_optimum_claim=True,no_independent_Pu_cohort=True,
       scripts_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (dest/'campaign_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
    if '--require-complete' in sys.argv:assert complete,'Campaign still has missing evaluation or figures'

if __name__=='__main__':run()
