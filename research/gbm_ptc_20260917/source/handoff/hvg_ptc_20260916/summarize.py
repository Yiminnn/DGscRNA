#!/usr/bin/env python3
"""Patient-aware summaries, paired effects and selection sensitivity from saved results."""
import argparse
import json
from pathlib import Path
from common import OUT, FEATURES, require_slurm, samples, sha, utc, write_json

METRICS=['partition_ari','partition_pair_f1','terminal_strict_L1_macroF1_present',
         'terminal_strict_L1_macroF1_fixed11','terminal_lfine_macroF1',
         'terminal_lfine_fixed22_setF1','terminal_coverage','terminal_strict_L1_accuracy','noise_rate']


def paired_delta(frame,feature,metric,cohort):
    import numpy as np
    from scipy.stats import wilcoxon
    keys=['sample','patient']
    a=frame[frame.feature.eq('all')][keys+[metric]].rename(columns={metric:'all_value'})
    b=frame[frame.feature.eq(feature)][keys+[metric]].rename(columns={metric:'hvg_value'})
    joined=a.merge(b,on=keys,validate='one_to_one')
    joined['delta']=joined.hvg_value-joined.all_value
    valid=joined.dropna(subset=['delta'])
    patient=valid.groupby('patient').delta.mean()
    result=dict(cohort=cohort,feature=feature,metric=metric,n_expected_samples=len(joined),
                n_paired_samples=len(valid),n_paired_patients=len(patient),
                n_expected_patients=joined.patient.nunique(),delta_mean=None,ci_lower=None,ci_upper=None,p_value=None)
    if len(patient):
        x=patient.to_numpy()
        rng=np.random.default_rng(20260916)
        boot=x[rng.integers(0,len(x),size=(10000,len(x)))].mean(1)
        result.update(delta_mean=float(x.mean()),delta_median=float(np.median(x)),
                      ci_lower=float(np.quantile(boot,.025)),ci_upper=float(np.quantile(boot,.975)),
                      n_patients_positive=int((x>0).sum()),n_patients_negative=int((x<0).sum()),
                      p_value=float(wilcoxon(x,zero_method='zsplit',method='auto').pvalue) if np.any(x) else 1.)
    valid=valid.copy();valid['feature']=feature;valid['metric']=metric;valid['cohort']=cohort
    return result,valid


def run(allow_partial=False):
    require_slurm()
    import numpy as np
    import pandas as pd
    from statsmodels.stats.multitest import multipletests
    out=OUT/'summary';out.mkdir(parents=True,exist_ok=True)
    frames=[];sourcehash={};incomplete=[]
    for sample in samples():
        p=OUT/'evaluation'/sample
        if not (p/'COMPLETE').exists():
            incomplete.append(sample)
            if allow_partial and (p/'metrics.csv').exists():frames.append(pd.read_csv(p/'metrics.csv'))
            continue
        m=json.loads((p/'manifest.json').read_text())
        assert m['status']=='completed' and m['expected']==330
        assert sha(p/'manifest.json')==(p/'COMPLETE').read_text().strip()
        assert sha(p/'metrics.csv')==m['outputs']['metrics.csv']
        sourcehash[sample]=sha(p/'manifest.json')
        frames.append(pd.read_csv(p/'metrics.csv'))
    if incomplete and not allow_partial:
        raise RuntimeError(f'Evaluation not complete for {len(incomplete)} samples: {incomplete}')
    if not frames:raise RuntimeError('No evaluated results')
    data=pd.concat(frames,ignore_index=True)
    assert not data.duplicated(['sample','condition']).any()
    if not incomplete:assert len(data)==39930 and data['sample'].nunique()==121
    for metric in METRICS:
        if metric not in data:data[metric]=np.nan
    data.to_csv(out/'all_sample_conditions.csv.gz',index=False)
    data.groupby(['feature','dr','clusterer','status'],dropna=False).size().rename('n').reset_index().to_csv(out/'availability.csv',index=False)
    cohorts={'primary97':data[data.evaluable.eq(True)],'all121':data,'excluded24':data[data.evaluable.eq(False)]}
    summaries=[]
    descriptors=['condition','geometry_id','arm_id','feature','dr','dim','input_space','seed','neighbors','min_dist',
                 'clusterer','k','covariance','min_cluster_size','min_samples','scoring_features','marker','cutoff','families']
    params=data[descriptors].drop_duplicates('condition').set_index('condition')
    for cname,d in cohorts.items():
        for cid,group in d.groupby('condition'):
            patient=group.groupby('patient')[METRICS].mean()
            row=dict(cohort=cname,condition=cid,n_expected_samples=len(group),n_expected_patients=group.patient.nunique(),
                     n_terminal_outputs=int(group.status.eq('completed').sum()),
                     n_structural_unavailable=int(group.status.str.startswith('structural').sum()),
                     n_numerical_unavailable=int(group.status.str.startswith('numerical').sum()),
                     n_training_executed=int(group.get('training_executed',False).eq(True).sum()),
                     **params.loc[cid].to_dict())
            for metric in METRICS:
                lo=-1 if metric=='partition_ari' else 0
                lower=group[['patient',metric]].copy();upper=lower.copy()
                lower[metric]=lower[metric].fillna(lo);upper[metric]=upper[metric].fillna(1)
                row.update({f'{metric}_patient_mean':patient[metric].mean(),
                            f'{metric}_sample_mean':group[metric].mean(),
                            f'{metric}_n_patients':int(patient[metric].notna().sum()),
                            f'{metric}_n_samples':int(group[metric].notna().sum()),
                            f'{metric}_fixed_cohort_lower':lower.groupby('patient')[metric].mean().mean(),
                            f'{metric}_fixed_cohort_upper':upper.groupby('patient')[metric].mean().mean()})
            summaries.append(row)
    pd.DataFrame(summaries).to_csv(out/'condition_summary.csv',index=False)
    # Same-sample paired effects; patients, not cells or random seeds, are independent units.
    tests=[];deltas=[]
    endpoints=['partition_ari','terminal_strict_L1_macroF1_present','terminal_lfine_macroF1']
    paths=[('direct_UMAP2','UMAP',2,'genes'),('PCA30_UMAP2','UMAP',2,'pca30'),('PCA30','PCA',30,'genes')]
    for cname,d in cohorts.items():
        for path,dr,dim,space in paths:
            subset=d[d.dr.eq(dr)&d.dim.eq(dim)&d.input_space.eq(space)&d.seed.eq(42)&
                     d.neighbors.eq(15)&d.min_dist.eq(.1)&d.clusterer.eq('HDBSCAN')&
                     d.min_cluster_size.eq(15)&d.min_samples.eq(15)&d.scoring_features.eq('all')&d.feature.isin(FEATURES)]
            for metric in endpoints:
                for feature in FEATURES[1:]:
                    row,delta=paired_delta(subset,feature,metric,cname)
                    row['path']=path;delta['path']=path;tests.append(row);deltas.append(delta)
    test=pd.DataFrame(tests)
    for _,idx in test.groupby(['cohort','path','metric']).groups.items():
        valid=test.loc[idx].p_value.notna()
        used=test.loc[idx][valid].index
        if len(used):test.loc[used,'p_holm_5_HVG_levels']=multipletests(test.loc[used,'p_value'],method='holm')[1]
    test['primary_contrast']=test.cohort.eq('primary97')&test.path.eq('direct_UMAP2')&test.feature.eq('hvg2000')&test.metric.isin(endpoints[:2])
    primary_idx=test.index[test.primary_contrast & test.p_value.notna()]
    if len(primary_idx):
        test.loc[primary_idx,'p_holm_2_primary_endpoints']=multipletests(test.loc[primary_idx,'p_value'],method='holm')[1]
    test.to_csv(out/'paired_patient_effects.csv',index=False)
    pd.concat(deltas,ignore_index=True).to_csv(out/'paired_sample_deltas.csv.gz',index=False)
    # The internal K criterion is evaluated without accessing reference scores for selection.
    ktable=data[data.families.str.contains('E2_K_selection')].copy()
    if len(ktable):
        selected=[]
        for (sample,feature,cl),d in ktable.groupby(['sample','feature','clusterer']):
            candidates=d[d.internal_silhouette.notna()].sort_values(['internal_silhouette','k'],ascending=[False,True])
            if len(candidates):
                choice=candidates.iloc[0].to_dict();choice['selection_rule']='max silhouette, tie smaller K; no reference labels'
                selected.append(choice)
        pd.DataFrame(selected).to_csv(out/'label_free_selected_K.csv',index=False)
    # Sensitivity to choosing HVG counts in other patients; held-out patients do not choose their own setting.
    cv=[]
    d=cohorts['primary97']
    fixed=d[d.dr.eq('UMAP')&d.dim.eq(2)&d.input_space.eq('genes')&d.seed.eq(42)&d.neighbors.eq(15)&
            d.min_dist.eq(.1)&d.clusterer.eq('HDBSCAN')&d.min_cluster_size.eq(15)&d.min_samples.eq(15)&
            d.scoring_features.eq('all')&d.feature.isin(FEATURES)]
    for metric in endpoints:
        table=fixed.pivot_table(index='patient',columns='feature',values=metric,aggfunc='mean').reindex(columns=FEATURES)
        for heldout in table.index:
            train=table.drop(index=heldout)
            allowed=train.notna().mean().ge(.9)
            score=train.mean().where(allowed)
            if score.notna().any():
                feature=score.idxmax()
                cv.append(dict(patient=heldout,metric=metric,chosen_feature=feature,
                    heldout_score=table.loc[heldout,feature],heldout_all=table.loc[heldout,'all'],
                    heldout_hvg2000=table.loc[heldout,'hvg2000'],training_patients=len(train),
                    selection='other-patient mean; >=90% training availability; fixed feature-order tie break'))
    pd.DataFrame(cv).to_csv(out/'leave_one_patient_out_HVG_selection.csv',index=False)
    # Seeds are algorithm variability, never additional biological replicates.
    seed=data[data.families.str.contains('E2_seeds')]
    if len(seed):
        seed.groupby(['sample','patient','feature','dr','dim','input_space','clusterer'],dropna=False)[endpoints].agg(['mean','std','min','max','count']).to_csv(out/'seed_variability.csv')
    coverage=[];feature_inventory=[]
    for sample in samples():
        prep=OUT/'prepared'/sample
        pm=json.loads((prep/'manifest.json').read_text())
        assert (prep/'PREPARED').read_text().strip()==sha(prep/'manifest.json')
        for feature,info in pm['feature_sets'].items():
            feature_inventory.append(dict(sample=sample,feature=feature,n_cells=pm['n_cells'],
                n_filtered_genes=pm['n_genes'],actual_geometry_genes=info['n_features'],status='available',
                final_DL_genes=pm['n_genes']))
        for feature,info in pm.get('unavailable_features',{}).items():
            feature_inventory.append(dict(sample=sample,feature=feature,n_cells=pm['n_cells'],
                n_filtered_genes=pm['n_genes'],actual_geometry_genes=np.nan,status=info['status'],
                final_DL_genes=pm['n_genes']))
        p=OUT/'prepared'/sample/'marker_coverage.csv'
        if p.exists():
            c=pd.read_csv(p);c['sample']=sample
            c['detected_marker_retention']=c.n_in_geometry/c.n_detected.replace(0,np.nan)
            coverage.append(c)
    pd.concat(coverage,ignore_index=True).to_csv(out/'marker_retention.csv.gz',index=False)
    pd.DataFrame(feature_inventory).to_csv(out/'actual_feature_inventory.csv',index=False)
    quality=[];quality_checks=[]
    for sample in samples():
        folder=OUT/'quality'/sample;p=folder/'geometry_quality.csv'
        if not p.exists():continue
        q=pd.read_csv(p)
        if (folder/'COMPLETE').exists():
            qm=json.loads((folder/'manifest.json').read_text())
            assert sha(folder/'manifest.json')==(folder/'COMPLETE').read_text().strip()
            assert qm['status']=='completed' and qm['n_complete']==qm['n_expected']==165
            assert len(q)==q.geometry_id.nunique()==165
            for name,digest in qm['outputs'].items():assert sha(folder/name)==digest
            for metric in ['trustworthiness_common_allgenes','neighbor_overlap_common_allgenes','distance_rv_common_allgenes']:
                vals=q[metric].dropna();assert vals.between(-1e-6,1+1e-6).all(),(sample,metric)
            identity=q[q.feature.eq('all')&q.dr.eq('none')]
            assert len(identity)==5 and set(identity.seed)=={42,7,17,29,101}
            # The untransformed all-gene representation is the diagnostic reference itself.
            assert identity.trustworthiness_common_allgenes.ge(.999).all()
            assert identity.neighbor_overlap_common_allgenes.ge(.99).all()
            assert identity.distance_rv_common_allgenes.abs().le(1e-6).all()
            quality_checks.append(dict(sample=sample,complete_representations=165,
                reference_self_n_seed_records=len(identity),
                reference_self_trustworthiness=identity.trustworthiness_common_allgenes.min(),
                reference_self_neighbor_overlap=identity.neighbor_overlap_common_allgenes.min(),
                reference_self_distance_rv=identity.distance_rv_common_allgenes.abs().max(),
                manifest_sha256=sha(folder/'manifest.json'),checks_passed=True))
        quality.append(q)
    if quality:pd.concat(quality,ignore_index=True).to_csv(out/'geometry_quality.csv.gz',index=False)
    pd.DataFrame(quality_checks).to_csv(out/'geometry_integrity_checks.csv',index=False)
    incomplete_quality=[s for s in samples() if not (OUT/'quality'/s/'COMPLETE').exists()]
    if incomplete_quality and not allow_partial:
        raise RuntimeError(f'Geometry diagnostics incomplete for {len(incomplete_quality)} samples')
    manifest=dict(timestamp=utc(),status='complete' if not incomplete and not incomplete_quality else 'partial',n_conditions=len(data),
                  n_samples=data['sample'].nunique(),evaluation_sources=sourcehash,incomplete_samples=incomplete,
                  incomplete_quality_samples=incomplete_quality,
                  summary_source_sha256=sha(Path(__file__)),
                  metrics_note='Terminal DL/refinement only; marker-only columns are labelled ablations. Structural missingness remains NA with fixed-denominator bounds.',
                  inference_note='Paired sample deltas averaged within patient; 10000 patient bootstrap replicates; Holm across five HVG levels per endpoint/path/cohort. Exploratory reanalysis of previously inspected data.',
                  outputs={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name not in ['manifest.json','COMPLETE']})
    write_json(out/'manifest.json',manifest)
    if not incomplete and not incomplete_quality:(out/'COMPLETE').write_text(sha(out/'manifest.json')+'\n')
    print(json.dumps({k:v for k,v in manifest.items() if k not in ['outputs','evaluation_sources']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--allow-partial',action='store_true')
    run(p.parse_args().allow_partial)
