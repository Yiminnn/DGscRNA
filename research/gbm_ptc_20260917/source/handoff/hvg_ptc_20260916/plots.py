#!/usr/bin/env python3
"""Saved-result publication figures; no models or reference-informed fitting."""
import argparse
import json
from common import OUT,FEATURES,require_slurm,sha,utc,write_json


def run(partial=False):
    require_slurm()
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    root=OUT/'summary';out=OUT/'figures';out.mkdir(parents=True,exist_ok=True)
    meta=json.loads((root/'manifest.json').read_text())
    if not partial:assert meta['status']=='complete' and (root/'COMPLETE').exists()
    data=pd.read_csv(root/'all_sample_conditions.csv.gz')
    data=data[data.evaluable.eq(True)].copy()
    summary=pd.read_csv(root/'condition_summary.csv')
    summary=summary[summary.cohort.eq('primary97')]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
        'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
    colors=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9']
    feature_order=['hvg500','hvg1000','hvg2000','hvg3000','hvg5000','all']
    labels=['500','1,000','2,000','3,000','5,000','All']
    captions=[]
    def save(fig,name,caption):
        if partial:fig.suptitle('PRELIMINARY: incomplete cohort',color='#A23B3B',fontsize=13)
        fig.tight_layout()
        for ext in ['png','pdf','svg']:fig.savefig(out/f'{name}.{ext}',dpi=240,bbox_inches='tight')
        plt.close(fig)
        captions.append(dict(name=name,caption=caption,status=meta['status'],
            files={ext:sha(out/f'{name}.{ext}') for ext in ['png','pdf','svg']}))
    def mean_ci(d,metric):
        vals=d.groupby('patient')[metric].mean().dropna().to_numpy()
        if not len(vals):return np.nan,np.nan,np.nan,0
        rng=np.random.default_rng(20260916)
        ci=np.quantile(vals[rng.integers(len(vals),size=(10000,len(vals)))].mean(1),[.025,.975])
        return vals.mean(),ci[0],ci[1],len(vals)
    def base(d):
        return d[d.seed.eq(42)&d.neighbors.eq(15)&d.min_dist.eq(.1)&d.clusterer.eq('HDBSCAN')&
                 d.min_cluster_size.eq(15)&d.min_samples.eq(15)&d.scoring_features.eq('all')]
    fixed=base(data)
    paths=[('Direct UMAP2','UMAP',2,'genes'),('PCA30 → UMAP2','UMAP',2,'pca30'),('PCA30','PCA',30,'genes')]
    endpoints=[('partition_ari','Partition ARI'),('terminal_strict_L1_macroF1_present','Terminal strict-L1 macro-F1'),
               ('terminal_lfine_macroF1','Terminal legacy Lfine concordance'),('noise_rate','Noise fraction')]
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for ax,(metric,title) in zip(axes.ravel(),endpoints):
        for i,(name,dr,dim,space) in enumerate(paths):
            d=fixed[fixed.dr.eq(dr)&fixed.dim.eq(dim)&fixed.input_space.eq(space)]
            vals=np.asarray([mean_ci(d[d.feature.eq(f)],metric)[:3] for f in feature_order])
            offset=(i-1)*.06
            ax.errorbar(np.arange(6)+offset,vals[:,0],yerr=np.maximum(0,np.vstack([vals[:,0]-vals[:,1],vals[:,2]-vals[:,0]])),
                        label=name,color=colors[i],marker='o',lw=1.5,capsize=2,markersize=4)
        ax.set_xticks(range(6),labels);ax.set_title(title,loc='left',weight='bold',fontsize=11)
        ax.set_xlabel('Genes used for geometry (categorical)');ax.grid(axis='y',alpha=.2)
    axes[0,0].legend(frameon=False,fontsize=9)
    save(fig,'hvg_response_curves','Primary eligible cohort. Each patient contributes equally after within-patient sample averaging. Error bars are 95% patient-bootstrap intervals; missing annotation conditions remain unavailable. All three paths use HDBSCAN15/15 and full-gene scoring/DL. The legacy Lfine score is set-valued concordance, not strict fine-label accuracy.')
    primary=summary[summary.families.str.contains('E1_primary_table')]
    dr_order=['PCA','FA','ICA','Isomap','UMAP','TSNE','none']
    cols=[f'{dr} / {c}' for dr in dr_order for c in ['KMeans','GMM','HDBSCAN']]
    primary=primary.copy();primary['method']=primary.dr+' / '+primary.clusterer
    fig,axes=plt.subplots(2,1,figsize=(15,7.5))
    for ax,metric,title in zip(axes,['partition_ari','terminal_strict_L1_macroF1_present'],
                              ['Partition ARI','Terminal strict-L1 macro-F1 (available conditions)']):
        a=primary.pivot(index='feature',columns='method',values=metric+'_patient_mean').reindex(index=feature_order,columns=cols)
        im=ax.imshow(a.to_numpy(),aspect='auto',cmap='cividis',vmin=-.1 if metric=='partition_ari' else 0,vmax=1)
        for y in range(len(a)):
            for x in range(len(a.columns)):
                v=a.iloc[y,x]
                if np.isfinite(v):ax.text(x,y,f'{v:.2f}',ha='center',va='center',fontsize=7.4,color='white' if v<.45 else '#172D3B')
                else:ax.text(x,y,'NA',ha='center',va='center',fontsize=7.4,color='#606060')
        ax.set_yticks(range(6),labels);ax.set_ylabel('Geometry genes');ax.set_title(title,loc='left',weight='bold')
        ax.set_xticks(range(len(cols)),cols,rotation=50,ha='right',fontsize=8)
        fig.colorbar(im,ax=ax,pad=.01,fraction=.015)
    save(fig,'primary_factorial_heatmaps','Six feature levels × seven representations × three clusterers. All reductions in this primary table are 2D; the no-reduction condition (none) uses the selected feature matrix directly. KMeans/GMM use predeclared K23; GMM covariance is diagonal. Patient means of available endpoints; availability and fixed-cohort bounds accompany this plot and are necessary for comparison.')
    method_path=OUT/'method_comparisons/paired_method_contrasts.csv'
    if method_path.exists():
        comparisons=pd.read_csv(method_path)
        comparisons=comparisons[comparisons.policy.eq('primary_rule')]
        comparators=[c for c in cols if c!='UMAP / HDBSCAN']
        metric_labels=[('partition_ari','Partition ARI'),('terminal_strict_L1_macroF1_present','Terminal strict-L1 macro-F1')]
        fig,axes=plt.subplots(1,2,figsize=(14,9),sharey=True)
        for ax,(metric,title) in zip(axes,metric_labels):
            for i,feature in enumerate(['all','hvg2000']):
                tab=comparisons[comparisons.feature.eq(feature)&comparisons.metric.eq(metric)].set_index('comparator').reindex(comparators)
                y=np.arange(len(tab))+(i-.5)*.25
                ax.errorbar(tab.delta_mean,y,xerr=np.maximum(0,np.vstack([tab.delta_mean-tab.ci_lower,tab.ci_upper-tab.delta_mean])),
                            fmt='o',markersize=4,capsize=2,color=colors[i],label=feature)
            ax.axvline(0,color='gray',lw=.9);ax.grid(axis='x',alpha=.15)
            ax.set_title(title,loc='left',weight='bold');ax.set_xlabel('Paired patient Δ: UMAP/HDBSCAN − comparator')
            ax.set_yticks(range(len(comparators)),comparators);ax.set_ylim(len(comparators)-.5,-.5)
        axes[0].legend(frameon=False)
        save(fig,'paired_method_comparisons','Exploratory comparisons of the frozen default factorial. Positive differences favor UMAP2/HDBSCAN15/15. Each contrast uses the same samples for both methods, followed by patient averaging; error bars are 95% patient-bootstrap intervals. All and HVG2000 are separate curves. Paired sample/patient counts and Holm corrections over all feature-method comparisons are saved in method_comparisons/paired_method_contrasts.csv. Pair availability can vary across comparators; these are not a common-cohort league table.')
        interaction=pd.read_csv(OUT/'method_comparisons/feature_method_interactions.csv')
        interaction=interaction[interaction.policy.eq('primary_rule')&interaction.feature.eq('hvg2000')]
        fig,axes=plt.subplots(1,2,figsize=(14,8.5),sharey=True)
        for ax,(metric,title) in zip(axes,metric_labels):
            tab=interaction[interaction.metric.eq(metric)].set_index('comparator').reindex(comparators)
            ax.errorbar(tab.delta_mean,np.arange(len(tab)),xerr=np.maximum(0,np.vstack([tab.delta_mean-tab.ci_lower,tab.ci_upper-tab.delta_mean])),
                        fmt='o',markersize=4,capsize=2,color=colors[0])
            ax.axvline(0,color='gray',lw=.9);ax.grid(axis='x',alpha=.15)
            ax.set_title(title,loc='left',weight='bold');ax.set_xlabel('HVG effect in UMAP/HDBSCAN − HVG effect in comparator')
            ax.set_yticks(range(len(comparators)),comparators);ax.set_ylim(len(comparators)-.5,-.5)
        save(fig,'feature_method_interaction','Exploratory four-condition interaction: [HVG2000 − all] in UMAP2/HDBSCAN15/15 minus [HVG2000 − all] in each comparator. Each sample must have all four scores; sample contrasts are then averaged by patient. Positive values indicate a larger HVG benefit in UMAP/HDBSCAN. Error bars are 95% patient-bootstrap intervals. Complete paired counts and multiplicity adjustments are supplied in the corresponding CSV; unavailable final annotations remain missing.')
    fig,axes=plt.subplots(1,2,figsize=(12,4.3),sharex=True)
    for ax,metric in zip(axes,['partition_ari','terminal_strict_L1_macroF1_present']):
        for i,space in enumerate(['genes','pca30']):
            for j,f in enumerate(['all','hvg2000']):
                d=fixed[fixed.dr.eq('UMAP')&fixed.input_space.eq(space)&fixed.feature.eq(f)]
                vals=np.asarray([mean_ci(d[d.dim.eq(n)],metric)[:3] for n in [2,10,30]])
                ax.errorbar([2,10,30],vals[:,0],yerr=np.maximum(0,np.vstack([vals[:,0]-vals[:,1],vals[:,2]-vals[:,0]])),
                            color=colors[2*i+j],marker=['o','s'][j],label=f'{"All" if f=="all" else "HVG2000"}; {"direct" if space=="genes" else "PCA30"}',capsize=2)
        ax.set_xticks([2,10,30]);ax.set_xlabel('UMAP dimensions used for clustering');ax.set_ylabel('Partition ARI' if metric=='partition_ari' else 'Terminal strict-L1 macro-F1');ax.grid(alpha=.15)
    axes[0].legend(frameon=False,fontsize=8)
    save(fig,'umap_dimension_factorial','Separate factorial sensitivity to all vs HVG2000, direct genes vs PCA30 input, and UMAP output dimension. All use HDBSCAN15/15. Display dimensionality is not substituted for the actual clustering input.')
    fig,axes=plt.subplots(1,2,figsize=(12,4.4))
    retention=pd.read_csv(root/'marker_retention.csv.gz')
    eligible=set(data['sample'].unique())
    r=retention[retention['sample'].isin(eligible)&retention.feature.isin(feature_order)]
    r=r.groupby(['sample','feature']).detected_marker_retention.mean().reset_index()
    for i,f in enumerate(feature_order):
        v=r.loc[r.feature.eq(f),'detected_marker_retention'].dropna()
        if len(v):axes[0].boxplot([v],positions=[i],widths=.55,showfliers=False,patch_artist=True,
                                  boxprops={'facecolor':'#D5E9F5','edgecolor':'#0072B2'},medianprops={'color':'#172D3B'})
    axes[0].set_xticks(range(6),labels);axes[0].set(ylabel='Mean per-panel detected-marker retention',xlabel='Geometry genes',ylim=(0,1.03))
    base2000=fixed[fixed.feature.eq('hvg2000')]
    pairs=[]
    metric='terminal_strict_L1_macroF1_present'
    for name,dr,dim,space in paths:
        b=base2000[base2000.dr.eq(dr)&base2000.dim.eq(dim)&base2000.input_space.eq(space)]
        c=data[data.feature.eq('hvg2000')&data.dr.eq(dr)&data.dim.eq(dim)&data.input_space.eq(space)&
               data.seed.eq(42)&data.neighbors.eq(15)&data.min_dist.eq(.1)&data.clusterer.eq('HDBSCAN')&
               data.min_cluster_size.eq(15)&data.min_samples.eq(15)&data.scoring_features.eq('geometry')]
        x=b[['sample','patient',metric]].merge(c[['sample',metric]],on='sample',suffixes=('_allscore','_hvgscore'),validate='one_to_one')
        x['delta']=x[metric+'_hvgscore']-x[metric+'_allscore'];x['path']=name;pairs.append(x)
    d=pd.concat(pairs,ignore_index=True);d.to_csv(out/'scoring_truncation_paired.csv',index=False)
    for i,(name,*_) in enumerate(paths):
        x=d[d.path.eq(name)].groupby('patient').delta.mean().dropna()
        axes[1].scatter(np.full(len(x),i)+np.random.default_rng(42).uniform(-.16,.16,len(x)),x,alpha=.55,s=12,color=colors[i])
        if len(x):axes[1].plot([i-.2,i+.2],[x.mean()]*2,color='black',lw=2)
    axes[1].axhline(0,color='gray',lw=.8);axes[1].set_xticks(range(3),[p[0] for p in paths],rotation=12)
    axes[1].set_ylabel('Patient Δ terminal F1: HVG scoring − full scoring')
    save(fig,'marker_retention_and_scoring_ablation','Left: marker retention in geometry, averaged equally over detected panels within each sample; primary scoring still sees full genes. Right: explicit scoring-truncation ablation at identical HVG2000 geometry. Points are patient means of paired sample differences; black bars are patient means. This distinguishes geometry selection from discarding scoring markers.')
    qpath=root/'geometry_quality.csv.gz'
    if qpath.exists():
        q=pd.read_csv(qpath).merge(data[['sample','patient']].drop_duplicates(),on='sample',how='inner',validate='many_to_one')
        q=q[q.dr.eq('UMAP')&q.seed.eq(42)&q.neighbors.eq(15)&q.min_dist.eq(.1)&q.feature.isin(['all','hvg2000'])]
        fig,axes=plt.subplots(1,3,figsize=(14,4))
        for ax,metric in zip(axes,['trustworthiness_common_allgenes','neighbor_overlap_common_allgenes','distance_rv_common_allgenes']):
            for i,(space,f) in enumerate((s,f) for s in ['genes','pca30'] for f in ['all','hvg2000']):
                d=q[q.input_space.eq(space)&q.feature.eq(f)]
                vals=np.asarray([mean_ci(d[d.dim.eq(n)],metric)[:3] for n in [2,10,30]])
                ax.plot([2,10,30],vals[:,0],marker='o',color=colors[i],label=f'{"All" if f=="all" else "HVG2000"}; {"direct" if space=="genes" else "PCA30"}')
            title={'trustworthiness_common_allgenes':'Trustworthiness','neighbor_overlap_common_allgenes':'Neighbor overlap','distance_rv_common_allgenes':'Distance residual variance'}[metric]
            ax.set_xticks([2,10,30]);ax.set_xlabel('UMAP dimensions');ax.set_title(title,fontsize=10)
            ax.grid(alpha=.15)
        axes[0].legend(frameon=False,fontsize=8)
        save(fig,'common_reference_geometry_diagnostics','All feature arms use the same all-gene scaled reference. Trustworthiness and neighborhood overlap use the same label-free random subset of at most 2000 cells, k=20. Distance residual variance uses the same 20000 sampled cell pairs; lower residual variance indicates better pair-distance preservation. Patient means shown.')
    secondary=OUT/'singleton_robustness_summary/secondary_primary_factorial.csv'
    if secondary.exists():
        second=pd.read_csv(secondary)
        second['method']=second.dr+' / '+second.clusterer
        fig,axes=plt.subplots(2,1,figsize=(15,7.3))
        for ax,field,title,maximum in zip(axes,['patient_mean','n_available_samples'],
            ['Secondary singleton-abstention rule: terminal strict-L1 macro-F1','Available sample count (primary cohort)'],[1,97]):
            tab=second.pivot(index='feature',columns='method',values=field).reindex(index=feature_order,columns=cols)
            im=ax.imshow(tab.to_numpy(),aspect='auto',cmap='cividis',vmin=0,vmax=maximum)
            for y in range(len(tab)):
                for x in range(len(tab.columns)):
                    v=tab.iloc[y,x]
                    if np.isfinite(v):ax.text(x,y,f'{v:.2f}' if field=='patient_mean' else str(int(v)),ha='center',va='center',fontsize=7.4,color='white' if v<maximum*.45 else '#172D3B')
            ax.set_yticks(range(6),labels);ax.set_xticks(range(len(cols)),cols,rotation=50,ha='right',fontsize=8)
            ax.set_title(title,loc='left',weight='bold');fig.colorbar(im,ax=ax,pad=.01,fraction=.015)
        save(fig,'secondary_singleton_robustness','Secondary implementation sensitivity, added after identifying singleton-DEG failures. All frozen cells, clusters, markers and DL settings remain fixed. Unsupported singleton groups seed Undecided; the same DL/refinement handles their cells. The primary-rule table is preserved separately. Remaining infeasible fixed DL splits stay unavailable. The lower panel displays sample availability, preventing a high mean from hiding failures.')
    if not partial:assert (OUT/'singleton_robustness_summary/COMPLETE').exists(), 'Secondary sensitivity verification pending'
    write_json(out/'figure_manifest.json',dict(timestamp=utc(),status=meta['status'],figures=captions,
        decision_tree='decision_tree_main.pdf',source_summary_sha256=sha(root/'manifest.json'),
        source_methods_sha256=sha(OUT/'method_comparisons/manifest.json') if (OUT/'method_comparisons/manifest.json').exists() else None,
        source_secondary_sha256=sha(OUT/'singleton_robustness_summary/manifest.json') if (OUT/'singleton_robustness_summary/manifest.json').exists() else None,
        plotting_source_sha256=sha(__file__),
        decision_tree_sha256={ext:sha(out/f'decision_tree_main.{ext}') for ext in ['pdf','svg','png']},
        tables={'scoring_truncation_paired.csv':sha(out/'scoring_truncation_paired.csv')}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--partial',action='store_true');run(p.parse_args().partial)
