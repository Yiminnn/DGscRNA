"""Summarize tested workflow choices and paired marker-by-DL module effects."""
import os
from common import OUT,require_slurm,checked,sha,write_json,complete,utc

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    from scipy.stats import wilcoxon
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    assert checked(OUT/'summary','aggregate_manifest.json','AGGREGATE_COMPLETE')
    assert checked(OUT/'controls_summary')
    dest=OUT/'workflow_choice_summary';dest.mkdir(exist_ok=True)
    if checked(dest):return
    source=OUT/'summary/marker_DL_factorial_heldout.csv'
    data=pd.read_csv(source)
    assert set(data.budget)=={'hvg2000'} and set(data.route)=={'UMAP2_HDBSCAN_R'}
    assert set(data.cutoff)=={'mean'} and set(data.seed)=={42}
    assert not data.duplicated(['cohort','patient','marker_selection','stage']).any()
    columns=[('fixed_marker','marker_only'),('fixed_marker','terminal090'),
             ('training_patient_selected_marker','marker_only'),('training_patient_selected_marker','terminal090')]
    # These are the four simple effects plus their difference-in-differences.
    # The descriptive outcomes already exist: this is retrospective inference,
    # not a prospective registration or an additional configuration search.
    definitions=[('DL_fixed_marker','DL with fixed marker',[-1,1,0,0]),
                 ('DL_selected_marker','DL with selected marker',[0,0,-1,1]),
                 ('marker_selection_initial','Marker selection before DL',[-1,0,1,0]),
                 ('marker_selection_terminal','Marker selection after DL',[0,-1,0,1]),
                 ('interaction','Marker-selection by DL interaction',[1,-1,-1,1])]
    rng=np.random.default_rng(20260917);rows=[];differences=[];point_rows=[]
    for cohort,g in data.groupby('cohort',sort=True):
        matrices={m:g.pivot(index='patient',columns=['marker_selection','stage'],values=m).reindex(columns=columns)
                  for m in ['macroF1_present','coverage']}
        patients=matrices['macroF1_present'].index
        assert len(patients)==(55 if cohort=='primary97' else 59)
        assert all(np.isfinite(t.to_numpy()).all() for t in matrices.values())
        draws=rng.integers(0,len(patients),size=(10000,len(patients)))
        for i,(marker,stage) in enumerate(columns):
            point_rows.append(dict(cohort=cohort,marker_selection=marker,stage=stage,n_patients=len(patients),
                macroF1_present_mean=float(matrices['macroF1_present'].iloc[:,i].mean()),
                coverage_mean=float(matrices['coverage'].iloc[:,i].mean())))
        for name,label,weights in definitions:
            record=dict(cohort=cohort,contrast=name,label=label,n_patients=len(patients))
            deltas={}
            for metric,table in matrices.items():
                z=table.to_numpy()
                simple={'DL_fixed_marker':z[:,1]-z[:,0], 'DL_selected_marker':z[:,3]-z[:,2],
                        'marker_selection_initial':z[:,2]-z[:,0], 'marker_selection_terminal':z[:,3]-z[:,1]}
                simple['interaction']=simple['DL_selected_marker']-simple['DL_fixed_marker']
                d=simple[name];deltas[metric]=d
                boot=d[draws].mean(axis=1)
                record[metric+'_mean_delta']=float(d.mean())
                record[metric+'_CI95_low']=float(np.quantile(boot,.025))
                record[metric+'_CI95_high']=float(np.quantile(boot,.975))
            d=deltas['macroF1_present']
            record['p_wilcoxon']=float(wilcoxon(d).pvalue) if np.any(np.abs(d)>1e-14) else 1.
            rows.append(record)
            for i,patient in enumerate(patients):
                differences.append(dict(cohort=cohort,contrast=name,patient=patient,
                    macroF1_present_delta=float(deltas['macroF1_present'][i]),coverage_delta=float(deltas['coverage'][i])))
    table=pd.DataFrame(rows);table['p_Holm']=1.
    for cohort,g in table.groupby('cohort'):
        ix=g.sort_values('p_wilcoxon',kind='stable').index
        table.loc[ix,'p_Holm']=np.minimum(1,np.maximum.accumulate(table.loc[ix,'p_wilcoxon'].to_numpy()*np.arange(5,0,-1)))
    table.to_csv(dest/'marker_DL_patient_paired.csv',index=False)
    pd.DataFrame(differences).to_csv(dest/'marker_DL_patient_differences.csv',index=False)
    points=pd.DataFrame(point_rows);points.to_csv(dest/'marker_DL_fourway_means.csv',index=False)
    old=pd.read_csv(OUT/'summary/marker_DL_factorial_summary.csv')
    parity=points.merge(old,on=['cohort','marker_selection','stage'],suffixes=('_new','_old'),validate='one_to_one')
    assert len(parity)==8
    for metric in ['macroF1_present_mean','coverage_mean']:
        np.testing.assert_allclose(parity[metric+'_new'],parity[metric+'_old'],rtol=0,atol=1e-14)
    selected=table[table.cohort=='primary97'].set_index('contrast').loc[[v[0] for v in definitions]]
    fig,axs=plt.subplots(1,2,figsize=(12,5.4),sharey=True,layout='constrained')
    y=np.arange(len(selected))
    for ax,metric,title in zip(axs,['macroF1_present','coverage'],['All-cell annotation macro-F1','Called-cell coverage']):
        means=selected[metric+'_mean_delta'].to_numpy()
        lo=selected[metric+'_CI95_low'].to_numpy();hi=selected[metric+'_CI95_high'].to_numpy()
        ax.hlines(y,lo,hi,color='#0072B2',lw=2);ax.scatter(means,y,color='#0072B2',s=36,zorder=3)
        ax.axvline(0,color='#888',lw=1);ax.set(title=title,xlabel='Paired mean change; 95% patient-bootstrap interval')
        ax.set_yticks(y,[v[1] for v in definitions]);ax.grid(axis='x',alpha=.15)
    axs[0].invert_yaxis()
    fig.suptitle('Marker choice and terminal DL have distinct effects\nOriginal HVG2000 / UMAP-HDBSCAN / mean cutoff; 55 heldout patients',fontsize=12)
    for ext in ['png','pdf']:fig.savefig(dest/f'marker_DL_patient_contrasts.{ext}',dpi=200,bbox_inches='tight')
    plt.close(fig)
    evidence=[]
    def add(node,comparison,row,scope,source_file,mean='mean_delta',low='CI95_low',high='CI95_high'):
        evidence.append(dict(node=node,comparison=comparison,n_patients=row.get('n_patients'),
            macroF1_mean_delta=row.get(mean),CI95_low=row.get(low),CI95_high=row.get(high),
            p_Holm=row.get('p_Holm'),scope_and_claim_limit=scope,source_file=source_file))
    fixed_path=OUT/'summary/fixed_marker_patient_paired.csv'
    fixed=pd.read_csv(fixed_path)
    fixed=fixed[(fixed.cohort=='primary97')&(fixed.route=='UMAP2_HDBSCAN_R')&(fixed.budget!='hvg2000')]
    for _,r in fixed.iterrows():
        add('Native feature budget',r.budget+' minus HVG2000',r,
            'Fixed historical glioma/mean on UMAP-HDBSCAN; geometry and DL genes change together, RNA scoring stays fixed. This does not prove that all genes are universally inferior.',
            'summary/fixed_marker_patient_paired.csv')
    geo_path=OUT/'controls_summary/geometry_vs_DL_paired_effects.csv';geo=pd.read_csv(geo_path)
    geo=geo[(geo.cohort=='primary97')&(geo.library=='CM2_glioma_other')&(geo.route=='UMAP2_HDBSCAN_R')]
    for _,r in geo.iterrows():
        if r.effect=='geometry_vs_2000_fixed_DL':
            node='Geometry genes alone';scope='Same full-RNA scoring and HVG2000 DL expression; only geometry feature budget changes.'
        else:
            node='DL input genes';scope='Same candidate geometry and scoring; candidate-budget DL genes minus fixed2000 DL genes.'
        add(node,r.budget+' / '+r.effect,r,scope,'controls_summary/geometry_vs_DL_paired_effects.csv')
    route_path=OUT/'summary/patient_heldout_route_comparisons.csv';routes=pd.read_csv(route_path)
    routes=routes[(routes.cohort=='primary97')&(routes.evidence=='database_or_external')]
    for _,r in routes.iterrows():
        add('Reduction and clustering route',r.route+' minus UMAP2_HDBSCAN_R',r,
            'Each route selects HVG/marker/cutoff on training patients. Conditional comparison of selected workflows, not a pure clusterer effect.',
            'summary/patient_heldout_route_comparisons.csv')
    for _,r in selected.reset_index().iterrows():
        add('Marker and DL modules',r.contrast,r,
            'Same mean cutoff, original geometry, patient folds and saved predictions. Selected marker chosen by training terminal scores for both stages. Marker-only remains an ablation.',
            'workflow_choice_summary/marker_DL_patient_paired.csv',
            'macroF1_present_mean_delta','macroF1_present_CI95_low','macroF1_present_CI95_high')
    for node,scope,source_file in [
        ('Normalization and input QC','Held fixed to the audited original-reference procedure; optimality versus other normalization/QC rules was not tested.','protocol/'),
        ('Density-scoring formula','Original density formula retained. Nonexistent alpha/epsilon mechanisms are not invented as ablations. No superiority over all alternative scoring formulas is established.','GBM_full_summary/METHODS_CORRECTIONS.md'),
        ('Embedding and clustering parameters','Limited sensitivity on three count-selected size pilots. Seeds are algorithm repeats, not patients; all-cohort optimal dimensions/minPts/resolution are not established.','controls_summary/representation_all_metrics.csv'),
        ('MLP architecture, epochs and seeds','54 controls on three count-selected pilots, with saved histories and separate initialization seeds. No claim of globally optimal architecture or epoch count.','controls_summary/MLP_all_metrics.csv')]:
        add(node,'Scope limit',{},scope,source_file)
    pd.DataFrame(evidence).to_csv(dest/'workflow_node_evidence.csv',index=False)
    note='''# Workflow-choice evidence and module contributions

All analysis reuses completed native-R predictions. No new annotation model is
fitted. The 2x2 table holds mean cutoff, HVG2000 and UMAP-HDBSCAN fixed; its selected
library is chosen on training-patient terminal scores and retained at both stages.
This marker-only selection is distinct from the39marker/cutoff choices used by
the comparator and Unknown analyses. Marker-only output is an ablation, never the
reported final DG-scRNA method output.

Five paired contrasts (two DL effects, two marker-selection effects, interaction)
are reported within each cohort. All patient differences remain available.
10000 patient resamples yield percentile intervals for mean effects; one common
resampling index per cohort keeps identical contrasts numerically identical.
Wilcoxon tests concern the signed-rank functional, not the mean. Holm correction
covers five macro-F1 tests per cohort; coverage is descriptive with mean intervals.
These are retrospective analyses of already visible descriptive results. Bootstrap
intervals condition on frozen cross-validated choices and do not refit selection.
The intervals and signed-rank tests need not give identical conclusions.

The workflow-node table carries through the original paired estimates and their
own correction families; it does not pool them into a new omnibus test or pretend
that every node was ablated. Zero marker-selection effects reflect identical
chosen outputs under the fixed cutoff, not proof that marker choice never matters.
Higher coverage alone is not evidence of better annotation. Uncertainty crossing
zero does not establish equivalence or noninferiority. Stability on three pilots
cannot establish optimal parameters across the full cohort. PTC integration,
marker retention and its historical group-specific branches are assessed separately.
'''
    (dest/'INTERPRETATION.md').write_text(note)
    write_json(dest/'manifest.json',dict(status='completed',no_new_fits=True,retrospective=True,
        n_paired_contrasts=len(table),n_node_evidence_rows=len(evidence),summary_mean_parity_passed=True,
        uncertainty='Paired patient bootstrap for mean effects; separate signed-rank tests; Holm five contrasts per cohort',
        source_files={str(p.relative_to(OUT)):sha(p) for p in [source,fixed_path,geo_path,route_path]},
        files={p.name:sha(p) for p in dest.iterdir() if p.suffix in ['.csv','.png','.pdf','.md']},
        job=os.environ['SLURM_JOB_ID'],source_sha256=sha(__file__),completed_at=utc()))
    complete(dest)

if __name__=='__main__':run()
