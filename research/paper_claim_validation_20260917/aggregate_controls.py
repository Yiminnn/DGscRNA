"""Separate geometry, DL-input, MLP and representation effects without re-fitting."""
import json
import os
from common import OUT, ROUTES, require_slurm, checked, complete, sha, write_json, utc

LIBRARIES = ['CM2_glioma_other', 'CM2_primary_all_context']

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    from scipy.stats import wilcoxon
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from dl_controls import CONFIGS
    from representation_controls import configurations
    import learning_curves

    assert checked(OUT/'summary', 'aggregate_manifest.json', 'AGGREGATE_COMPLETE')
    dest = OUT/'controls_summary'; dest.mkdir(exist_ok=True)
    cohort = pd.read_csv(OUT/'protocol/cohort.csv')
    pilots = (OUT/'protocol/pilot_samples.txt').read_text().split()
    sources = {}; geometry = []; mlp = []; rep = []; figures = []
    def read(path):
        assert checked(path), str(path)
        sources[str(path.relative_to(OUT))] = sha(path/'manifest.json')
        return pd.read_csv(path/'metrics.csv', dtype={'cutoff': str})
    for sample in cohort['sample']:
        for budget in ['hvg2000', 'hvg5000', 'all']:
            geometry.append(read(OUT/'GBM'/sample/budget/'evaluation_geometry_only_DL2000'))
    for sample in pilots:
        for budget in ['hvg2000', 'hvg5000']:
            for config in CONFIGS:
                mlp.append(read(OUT/'GBM_DL_controls'/sample/budget/config))
            for cfg in configurations(sample, budget):
                path = OUT/'GBM_representation_controls'/sample/budget/cfg['name']
                rep.append(read(path))
                assert (path/'clustering_and_terminal.png').exists() and (path/'clustering_and_terminal.pdf').exists()
                figures.append(dict(**cfg, png=str(path.relative_to(OUT)/'clustering_and_terminal.png'),
                                    pdf=str(path.relative_to(OUT)/'clustering_and_terminal.pdf')))
    assert (len(geometry), len(mlp), len(rep)) == (363, 54, 180)
    geo = pd.concat(geometry, ignore_index=True)
    m = pd.concat(mlp, ignore_index=True)
    r = pd.concat(rep, ignore_index=True)
    native = pd.read_csv(OUT/'summary/all_annotation_metrics.csv.gz', dtype={'cutoff': str})
    native = native[(native.library.isin(LIBRARIES)) & (native.cutoff=='mean') &
                    native.budget.isin(['hvg2000','hvg5000','all'])].copy()
    measures = ['macroF1_present','accuracy','coverage','unknown_rate']
    # At 2000 the native and geometry-only experiments must coincide exactly.
    keys = ['sample','budget','route','library','cutoff','stage']
    same = geo[geo.budget=='hvg2000'].merge(native[native.budget=='hvg2000'], on=keys,
                                             suffixes=('_geometry','_native'), validate='one_to_one')
    assert len(same)==121*4*2*3
    for metric in measures:
        assert np.array_equal(same[metric+'_geometry'].to_numpy(), same[metric+'_native'].to_numpy()), metric
    joint = pd.concat([native, geo], ignore_index=True)
    geo.to_csv(dest/'geometry_only_all_metrics.csv.gz', index=False, compression='gzip')
    m.to_csv(dest/'MLP_all_metrics.csv', index=False)
    r.to_csv(dest/'representation_all_metrics.csv', index=False)
    pd.DataFrame(figures).to_csv(dest/'representation_figure_index.csv', index=False)
    spec = ['family','budget','route','library','stage']
    patient_frames = []
    for name, frame in [('primary97',joint[joint.primary]), ('all121',joint)]:
        patient = frame.groupby(['patient']+spec)[measures].mean().reset_index()
        patient.insert(0,'cohort',name); patient_frames.append(patient)
    patient = pd.concat(patient_frames, ignore_index=True)
    patient.to_csv(dest/'geometry_and_native_patient.csv', index=False)
    avg = patient.groupby(['cohort']+spec)[measures].agg(['mean','std','count'])
    avg.columns = ['_'.join(c) for c in avg.columns]
    avg.reset_index().to_csv(dest/'geometry_and_native_summary.csv', index=False)
    final = patient[patient.stage=='terminal090']
    rng = np.random.default_rng(20260917); contrasts = []
    for (name, lib, route), group in final.groupby(['cohort','library','route']):
        g = group[group.family=='geometry_only_fixed_DL2000']
        base = g[g.budget=='hvg2000'][['patient']+measures]
        for budget in ['hvg5000','all']:
            candidate = g[g.budget==budget][['patient']+measures]
            whole = group[(group.family=='native_R_budget') & (group.budget==budget)][['patient']+measures]
            for effect, left, right in [('geometry_vs_2000_fixed_DL',candidate,base),
                                         ('DL_budget_at_same_geometry',whole,candidate)]:
                p = left.merge(right,on='patient',suffixes=('_candidate','_baseline'),validate='one_to_one')
                delta = (p.macroF1_present_candidate-p.macroF1_present_baseline).to_numpy()
                boot = delta[rng.integers(0,len(delta),size=(10000,len(delta)))].mean(axis=1)
                contrasts.append(dict(cohort=name,library=lib,route=route,budget=budget,effect=effect,
                    n_patients=len(delta),mean_delta=float(delta.mean()),CI95_low=float(np.quantile(boot,.025)),
                    CI95_high=float(np.quantile(boot,.975)),
                    p=float(wilcoxon(delta).pvalue) if np.any(np.abs(delta)>1e-14) else 1.))
    contrasts = pd.DataFrame(contrasts); contrasts['p_Holm']=1.
    for _, frame in contrasts.groupby(['cohort','effect']):
        ix=frame.sort_values('p').index
        contrasts.loc[ix,'p_Holm']=np.minimum(1,np.maximum.accumulate(contrasts.loc[ix,'p'].to_numpy()*np.arange(len(ix),0,-1)))
    contrasts.to_csv(dest/'geometry_vs_DL_paired_effects.csv',index=False)
    fig, axs = plt.subplots(2,2,figsize=(12,8),layout='constrained')
    colors = ['#0072B2','#D55E00','#009E73','#CC79A7']
    for i, lib in enumerate(LIBRARIES):
        for j, family in enumerate(['native_R_budget','geometry_only_fixed_DL2000']):
            ax=axs[i,j]
            view=final[(final.cohort=='primary97')&(final.library==lib)&(final.family==family)]
            for route,color in zip(ROUTES,colors):
                values=view[view.route==route].groupby('budget').macroF1_present.mean().reindex(['hvg2000','hvg5000','all'])
                ax.plot(range(3),values,marker='o',color=color,label=route)
            ax.set_xticks(range(3),['2,000','5,000','All']);ax.set_ylim(0,1)
            ax.set(title=lib+'\n'+('Geometry + DL feature budget' if j==0 else 'Geometry only; fixed DL2000'),
                   xlabel='Geometry feature budget',ylabel='Patient-weighted terminal macro-F1')
    axs[0,0].legend(fontsize=7,frameon=False)
    for ext in ['png','pdf']:fig.savefig(dest/f'geometry_and_DL.{ext}',dpi=210,bbox_inches='tight')
    plt.close(fig)
    # Algorithm repeats on the same three samples are not biological replicates.
    anchor=native[native['sample'].isin(pilots)&native.budget.isin(['hvg2000','hvg5000'])&native.stage.eq('terminal090')].copy()
    anchor['stage']='final090';anchor['space']=anchor.route.str.replace('_SNN','',regex=False).str.replace('_HDBSCAN_R','',regex=False)
    anchor['method']=np.where(anchor.route.str.contains('HDBSCAN'),'HDBSCAN_R','SNN')
    anchor['minPts']=50;anchor['resolution']=.5;anchor['embedding_seed']=42;anchor['kind']='anchor';anchor['name']='original_default'
    representations=pd.concat([r[r.stage=='final090'],anchor],ignore_index=True)
    representations.to_csv(dest/'representation_with_original_anchors.csv',index=False)
    for title, frame, xkey, filename in [
        ('MLP controls: clusters, marker seeds and split fixed',m[m.stage=='final090'],'control','MLP_controls'),
        ('Dimension / clustering parameter / embedding-seed controls',representations,'name','representation_controls')]:
        fig,axs=plt.subplots(3,2,figsize=(16,12),layout='constrained')
        order=list(CONFIGS) if xkey=='control' else list(dict.fromkeys(frame[xkey]))
        for i,sample in enumerate(pilots):
            for j,budget in enumerate(['hvg2000','hvg5000']):
                ax=axs[i,j];view=frame[(frame['sample']==sample)&(frame.budget==budget)]
                for k,lib in enumerate(LIBRARIES):
                    styles=list(zip(ROUTES,['o','^','s','D'])) if xkey=='control' else [('SNN','o'),('HDBSCAN_R','^')]
                    for method,marker in styles:
                        mask=view.route.eq(method) if xkey=='control' else view.route.str.endswith('_'+method)
                        take=view[view.library.eq(lib)&mask]
                        ax.scatter([order.index(v)+(k-.5)*.18 for v in take[xkey]],take.macroF1_present,
                                   s=15,marker=marker,alpha=.7,color=['#0072B2','#D55E00'][k],label=lib+' / '+method)
                ax.set_xticks(range(len(order)),order,rotation=90,fontsize=5 if xkey=='name' else 7)
                ax.set(title=f'{sample} / {budget}',ylabel='Terminal macro-F1',ylim=(0,1))
        axs[0,0].legend(fontsize=6,frameon=False,ncol=2)
        fig.suptitle(title+'\nPrespecified size pilots; points are algorithm conditions, not additional patients')
        for ext in ['png','pdf']:fig.savefig(dest/f'{filename}.{ext}',dpi=190,bbox_inches='tight')
        plt.close(fig)
    learning_curves.run(dest)
    write_json(dest/'manifest.json',dict(status='completed',geometry_units=363,MLP_units=54,representation_units=180,
        identical_2000_anchor=True,source_manifests=sources,
        learning_curves_manifest_sha256=sha(dest/'learning_curves_manifest.json'),
        scopes=dict(geometry='All121samples/59patients; primary97/55separate',
                    MLP_and_representation='Three count-selected samples, descriptive stability only; no all-cohort optimality inference'),
        files={p.name:sha(p) for p in dest.iterdir() if p.suffix in ['.csv','.gz','.png','.pdf']},
        job=os.environ['SLURM_JOB_ID'],source_sha256=sha(__file__),completed_at=utc()))
    complete(dest)

if __name__=='__main__':run()
