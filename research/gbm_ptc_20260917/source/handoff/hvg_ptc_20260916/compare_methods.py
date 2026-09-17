#!/usr/bin/env python3
"""Exploratory paired method and feature-by-method contrasts of frozen outputs."""
import json
from pathlib import Path
from common import OUT,FEATURES,require_slurm,sha,utc,write_json


def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    from scipy.stats import wilcoxon
    from statsmodels.stats.multitest import multipletests
    assert (OUT/'summary/COMPLETE').exists()
    assert (OUT/'singleton_robustness_summary/COMPLETE').exists()
    out=OUT/'method_comparisons';out.mkdir(exist_ok=True)
    sources={'primary_rule':OUT/'summary/all_sample_conditions.csv.gz',
             'secondary_singleton_rule':OUT/'singleton_robustness_summary/secondary_all_conditions.csv.gz'}
    tables={name:pd.read_csv(path) for name,path in sources.items()}
    primary=tables['primary_rule']
    assert len(primary)==39930
    tests=[];interactions=[];ranking=[];dl=[]
    endpoints=['partition_ari','terminal_strict_L1_macroF1_present']
    baseline='UMAP / HDBSCAN'
    def stats(delta,expected):
        delta=delta.dropna()
        patient=delta.groupby(level='patient').mean()
        r=dict(n_expected_samples=expected,n_paired_samples=len(delta),n_paired_patients=len(patient),
               delta_mean=np.nan,ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan)
        if len(patient):
            x=patient.to_numpy();rng=np.random.default_rng(20260916)
            boot=x[rng.integers(len(x),size=(10000,len(x)))].mean(1)
            r.update(delta_mean=float(x.mean()),ci_lower=float(np.quantile(boot,.025)),
                     ci_upper=float(np.quantile(boot,.975)),
                     p_value=float(wilcoxon(x,zero_method='zsplit',method='auto').pvalue) if np.any(x) else 1.)
        return r
    for policy,full in tables.items():
        d=full[full.evaluable.eq(True)&full.families.str.contains('E1_primary_table')].copy()
        assert d['sample'].nunique()==97 and len(d)==97*126
        d['method']=d.dr+' / '+d.clusterer
        for metric in endpoints:
            if policy=='secondary_singleton_rule' and metric=='partition_ari':continue
            wide=d.pivot(index=['sample','patient'],columns=['feature','method'],values=metric)
            methods=sorted(d.method.unique());assert len(methods)==21
            # Complete-case ranking uses the same samples/patients across all 21 methods.
            for feature in FEATURES:
                tab=wide[feature]
                shared=tab.dropna()
                shared_patients=shared.groupby(level='patient').mean()
                for method in methods:
                    available=tab[method].dropna().groupby(level='patient').mean()
                    ranking.append(dict(policy=policy,feature=feature,metric=metric,method=method,
                        n_available_samples=int(tab[method].notna().sum()),n_available_patients=len(available),
                        available_patient_mean=available.mean(),n_shared_samples=len(shared),
                        n_shared_patients=len(shared_patients),shared_patient_mean=shared_patients[method].mean(),
                        fixed_cohort_lower=tab[method].fillna(-1 if metric=='partition_ari' else 0).groupby(level='patient').mean().mean(),
                        fixed_cohort_upper=tab[method].fillna(1).groupby(level='patient').mean().mean()))
                    if method==baseline:continue
                    delta=tab[baseline]-tab[method]
                    tests.append(dict(policy=policy,feature=feature,metric=metric,comparator=method,
                        contrast='UMAP/HDBSCAN minus comparator',**stats(delta,len(tab))))
                    if feature=='all':continue
                    # Four-condition complete pairs isolate a feature-by-method interaction.
                    interaction=(wide[(feature,baseline)]-wide[('all',baseline)])-(wide[(feature,method)]-wide[('all',method)])
                    interactions.append(dict(policy=policy,feature=feature,metric=metric,comparator=method,
                        contrast='HVG effect in UMAP/HDBSCAN minus HVG effect in comparator',**stats(interaction,len(tab))))
    for rows,name in [(tests,'paired_method_contrasts.csv'),(interactions,'feature_method_interactions.csv')]:
        table=pd.DataFrame(rows)
        for _,index in table.groupby(['policy','metric']).groups.items():
            used=table.loc[index].dropna(subset=['p_value']).index
            table.loc[used,'p_holm_all_contrasts_in_policy_endpoint']=multipletests(table.loc[used,'p_value'],method='holm')[1]
        table.to_csv(out/name,index=False)
    pd.DataFrame(ranking).to_csv(out/'method_rankings_with_shared_cohort.csv',index=False)
    # Initial calls are explicitly an ablation. The terminal output remains the endpoint.
    d=primary[primary.evaluable.eq(True)&primary.families.str.contains('E1_primary_table')]
    for (feature,dr,clusterer),g in d.groupby(['feature','dr','clusterer']):
        index=g.set_index(['sample','patient'])
        delta=index.terminal_strict_L1_macroF1_present-index.marker_only_ablation_strict_L1_macroF1_present
        dl.append(dict(feature=feature,dr=dr,clusterer=clusterer,contrast='terminal minus marker-only ablation',
                       n_training_executed=int(g.training_executed.eq(True).sum()),**stats(delta,len(g))))
    pd.DataFrame(dl).to_csv(out/'terminal_refinement_effect.csv',index=False)
    # Predeclared A5 interventions, estimated from the unchanged saved conditions.
    d=primary[primary.evaluable.eq(True)&primary.seed.eq(42)&primary.neighbors.eq(15)&
        primary.min_dist.eq(.1)&primary.clusterer.eq('HDBSCAN')&primary.min_cluster_size.eq(15)&primary.min_samples.eq(15)]
    a5=[]
    for path,dr,dim,space in [('direct_UMAP2','UMAP',2,'genes'),('PCA30_UMAP2','UMAP',2,'pca30'),('PCA30','PCA',30,'genes')]:
        branch=d[d.dr.eq(dr)&d.dim.eq(dim)&d.input_space.eq(space)]
        reference=branch[branch.feature.eq('hvg2000')&branch.scoring_features.eq('all')].set_index(['sample','patient'])
        assert len(reference)==97 and reference.index.is_unique
        for intervention,feature,scoring in [('marker_union','hvg2000_markers','all'),('scoring_truncation','hvg2000','geometry')]:
            treatment=branch[branch.feature.eq(feature)&branch.scoring_features.eq(scoring)].set_index(['sample','patient'])
            assert len(treatment)==97 and set(treatment.index)==set(reference.index)
            for metric in endpoints:
                delta=treatment[metric]-reference[metric]
                if intervention=='scoring_truncation' and metric=='partition_ari':
                    assert delta.dropna().eq(0).all(), 'Scoring-only intervention changed partition scores'
                    continue
                a5.append(dict(intervention=intervention,path=path,metric=metric,
                    contrast='intervention minus HVG2000 geometry with full-gene scoring',**stats(delta,97)))
    a5=pd.DataFrame(a5)
    for _,index in a5.groupby(['intervention','metric']).groups.items():
        used=a5.loc[index].dropna(subset=['p_value']).index
        if len(used):a5.loc[used,'p_holm_3_paths']=multipletests(a5.loc[used,'p_value'],method='holm')[1]
    a5.to_csv(out/'marker_union_and_scoring_truncation_effects.csv',index=False)
    preprocessing=[]
    for dimension in [2,10,30]:
        branch=d[d.dr.eq('UMAP')&d.dim.eq(dimension)&d.feature.isin(['all','hvg2000'])&d.scoring_features.eq('all')]
        for metric in endpoints:
            wide=branch.pivot(index=['sample','patient'],columns=['feature','input_space'],values=metric)
            assert len(wide)==97 and len(wide.columns)==4
            delta=(wide[('hvg2000','genes')]-wide[('all','genes')])-(wide[('hvg2000','pca30')]-wide[('all','pca30')])
            preprocessing.append(dict(umap_dim=dimension,metric=metric,
                contrast='HVG2000 effect with direct genes minus HVG2000 effect after PCA30',**stats(delta,97)))
    preprocessing=pd.DataFrame(preprocessing)
    for _,index in preprocessing.groupby('metric').groups.items():
        used=preprocessing.loc[index].dropna(subset=['p_value']).index
        if len(used):preprocessing.loc[used,'p_holm_3_dimensions']=multipletests(preprocessing.loc[used,'p_value'],method='holm')[1]
    preprocessing.to_csv(out/'hvg_pca_interactions.csv',index=False)
    write_json(out/'manifest.json',dict(status='completed',timestamp=utc(),
        analysis_scope='Exploratory comparisons of the full frozen factorial; no refitting or outcome-guided parameter changes.',
        inference='Paired sample differences averaged by patient; 10000 patient bootstrap draws. Holm over all method contrasts within policy and endpoint. Four-condition complete pairs for interactions. Predeclared A5 marker/scoring interventions use Holm over three paths per intervention/endpoint. HVG-by-PCA interactions describe the predeclared UMAP factorial, with Holm over three dimensions per endpoint.',
        sources={k:sha(v) for k,v in sources.items()},source_sha256=sha(Path(__file__)),
        outputs={p.name:sha(p) for p in out.glob('*.csv')}))
    (out/'COMPLETE').write_text(sha(out/'manifest.json')+'\n')
    print('Verified exploratory method comparisons saved',flush=True)


if __name__=='__main__':run()
