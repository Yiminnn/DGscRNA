"""Exploratory patient-paired expression contrasts for terminal Unknown cells."""
import json
import os
from common import OUT,L1,require_slurm,checked,sha,write_json,complete,utc

MIN_CELLS=20
MIN_PATIENTS=6
BOOTSTRAPS=2000

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    from scipy.stats import binom
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from prediction_helpers import read_native,UNKNOWN
    assert checked(OUT/'comparison_summary')
    dest=OUT/'unknown_expression';dest.mkdir(exist_ok=True)
    cohort=pd.read_csv(OUT/'protocol/cohort.csv');cohort=cohort[cohort.primary]
    cv=pd.read_csv(OUT/'comparison_summary/patient_heldout_results.csv',dtype={'cutoff':str})
    cv=cv[(cv.method=='DG-scRNA')&(cv.cohort=='primary97')]
    genes=set();local={};sources={}
    for row in cohort.itertuples():
        p=OUT/'GBM'/row.sample/'hvg2000/Seurat_gene_names.csv'
        g=pd.read_csv(p,dtype=str,keep_default_na=False)
        assert g.Seurat.is_unique
        local[row.sample]=g;genes.update(g.Seurat)
    genes=sorted(genes);gene_index={g:i for i,g in enumerate(genes)};ng=len(genes)
    accumulated={};support=[];normalization=[]
    conditions=['fixed_glioma','training_patient_selected']
    for row in cohort.itertuples():
        src=OUT/'inputs'/row.sample;im=json.loads((src/'input_manifest.json').read_text())
        assert checked(src,'input_manifest.json','INPUT_COMPLETE')
        assert all(sha(src/n)==im['fitting_files'][n] for n in ['x.bin','i.bin','p.bin','genes.csv','cells_fit.csv'])
        X=sp.csr_matrix((np.fromfile(src/'x.bin',dtype='<f8'),np.fromfile(src/'i.bin',dtype='<i4'),
            np.fromfile(src/'p.bin',dtype='<i4')),shape=(im['n_cells'],im['n_genes']))
        X.eliminate_zeros();assert (X.data>0).all()
        g=local[row.sample]
        assert list(pd.read_csv(src/'genes.csv',dtype=str).gene)==list(g.source)
        truth=pd.read_csv(OUT/'evaluation_inputs'/row.sample/'truth.csv.gz',dtype=str,keep_default_na=False)
        fit=pd.read_csv(src/'cells_fit.csv',dtype=str,keep_default_na=False)
        assert list(fit.cell_id)==list(truth.cell_id)
        totals=np.asarray(X.sum(axis=1)).ravel();assert (totals>0).all()
        X.data=np.log1p(X.data*np.repeat(10000/totals,np.diff(X.indptr)))
        indices=np.asarray([gene_index[v] for v in g.Seurat])
        # Independently verify the normalization against the R-exported DL matrix.
        prep=OUT/'GBM'/row.sample/'hvg2000';features=(prep/'DL_features.txt').read_text().splitlines()
        pos={v:i for i,v in enumerate(g.Seurat)};selected=[pos[v] for v in features]
        dl=np.memmap(prep/'DL.float32.bin',dtype='<f4',mode='r',shape=(len(truth),len(features)))
        cells=np.unique(np.linspace(0,len(truth)-1,min(16,len(truth)),dtype=int))
        got=X[cells][:,selected].toarray().astype(np.float32)
        np.testing.assert_allclose(got,dl[cells],rtol=1e-6,atol=1e-6)
        normalization.append(dict(sample=row.sample,n_checked_cells=len(cells),n_checked_genes=len(features),
            max_absolute_difference=float(np.max(np.abs(got-dl[cells])))))
        del dl,got
        chosen=cv[cv.patient==row.patient].iloc[0]
        for condition,library,cutoff in [('fixed_glioma','CM2_glioma_other','mean'),
                                        ('training_patient_selected',chosen.library,chosen.cutoff)]:
            pred,path=read_native('DG-scRNA',row.sample,'hvg2000','UMAP2_HDBSCAN_R',library,cutoff)
            assert list(pred.cell_id)==list(truth.cell_id)
            sources[str(path.relative_to(OUT))]=sha(path)
            unknown=pred.prediction.isin(UNKNOWN).to_numpy()
            for label in L1:
                same=truth.L1.eq(label).to_numpy();a=same&unknown;b=same&~unknown
                na=int(a.sum());nb=int(b.sum());eligible=min(na,nb)>=MIN_CELLS
                support.append(dict(sample=row.sample,patient=row.patient,condition=condition,label=label,
                    n_Unknown=na,n_called=nb,eligible_pair=eligible,library=library,cutoff=cutoff))
                if not eligible:continue
                values=[]
                for mask in [a,b]:
                    block=X[mask]
                    values.extend([np.asarray(block.mean(axis=0)).ravel(),block.getnnz(axis=0)/block.shape[0]])
                key=(condition,label,row.patient)
                if key not in accumulated:accumulated[key]=(np.zeros((4,ng)),np.zeros(ng,dtype=np.int16))
                sums,counts=accumulated[key]
                sums[:,indices]+=np.asarray(values);counts[indices]+=1
        sources[str((src/'input_manifest.json').relative_to(OUT))]=sha(src/'input_manifest.json')
        print('UNKNOWN_EXPRESSION_SAMPLE',row.sample,flush=True)
    pd.DataFrame(support).to_csv(dest/'sample_class_pair_support.csv',index=False)
    pd.DataFrame(normalization).to_csv(dest/'normalization_R_parity.csv',index=False)
    rng=np.random.default_rng(20260917);frames=[];summaries=[]
    def bh(p):
        order=np.argsort(p,kind='stable');q=np.empty(len(p))
        q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
        return q
    for condition in conditions:
        for label in L1:
            keys=sorted(k for k in accumulated if k[:2]==(condition,label));patients=[k[2] for k in keys]
            if not keys:
                summaries.append(dict(condition=condition,label=label,n_patients_with_pairs=0,n_tested_genes=0,status='no_eligible_pairs'))
                continue
            arrays=[]
            for key in keys:
                sums,counts=accumulated[key]
                arrays.append(np.divide(sums,counts[None,:],out=np.full_like(sums,np.nan),where=counts[None,:]>0))
            arr=np.asarray(arrays);delta=arr[:,0,:]-arr[:,2,:]
            valid=np.isfinite(delta);n=valid.sum(0)
            avg=np.divide(np.nansum(arr,axis=0),n[None,:],out=np.full((4,ng),np.nan),where=n[None,:]>0)
            eligible=(n>=MIN_PATIENTS)&(np.maximum(avg[1],avg[3])>=.1)
            name=condition+'__'+label.replace(' ','_')
            np.savez_compressed(dest/(name+'_patient_expression.npz'),patients=np.asarray(patients),genes=np.asarray(genes),
                mean_log1p_Unknown=arr[:,0,:],fraction_detected_Unknown=arr[:,1,:],
                mean_log1p_called=arr[:,2,:],fraction_detected_called=arr[:,3,:],delta=delta)
            ix=np.where(eligible)[0]
            summaries.append(dict(condition=condition,label=label,n_patients_with_pairs=len(keys),n_tested_genes=len(ix),
                status='tested' if len(ix) else 'insufficient_gene_or_patient_support'))
            if not len(ix):continue
            d=delta[:,ix];positive=(d>1e-12).sum(0);negative=(d< -1e-12).sum(0);effective=positive+negative
            p=np.minimum(1,2*binom.cdf(np.minimum(positive,negative),effective,.5))
            p[effective==0]=1
            # Resample patients, preserving class/condition pairing and observed
            # gene availability; no cell is an inferential replicate.
            weights=rng.multinomial(len(keys),np.full(len(keys),1/len(keys)),size=BOOTSTRAPS)
            lo=np.zeros(len(ix));hi=lo.copy()
            for first in range(0,len(ix),256):
                block=d[:,first:first+256];ok=np.isfinite(block)
                denominator=weights@ok.astype(float)
                bootstrap=np.divide(weights@np.nan_to_num(block),denominator,
                    out=np.full(denominator.shape,np.nan),where=denominator>0)
                lo[first:first+256],hi[first:first+256]=np.nanquantile(bootstrap,[.025,.975],axis=0)
            f=pd.DataFrame(dict(condition=condition,label=label,gene=np.asarray(genes)[ix],n_paired_patients=n[ix],
                n_positive=positive,n_negative=negative,n_tied=n[ix]-effective,
                mean_log1p_Unknown=avg[0,ix],mean_log1p_called=avg[2,ix],
                mean_fraction_detected_Unknown=avg[1,ix],mean_fraction_detected_called=avg[3,ix],
                mean_log1p_difference=np.nanmean(d,axis=0),median_log1p_difference=np.nanmedian(d,axis=0),
                CI95_mean_low=lo,CI95_mean_high=hi,sign_test_p=p,q_BH_within_contrast=bh(p)))
            frames.append(f)
    summaries=pd.DataFrame(summaries);summaries.to_csv(dest/'class_condition_support.csv',index=False)
    assert frames,'No paired expression contrast met the frozen minimum support'
    table=pd.concat(frames,ignore_index=True);table['q_BH_all_tested_contrasts']=bh(table.sign_test_p.to_numpy())
    table.to_csv(dest/'all_patient_paired_gene_contrasts.csv.gz',index=False)
    ranked=table.assign(abs_effect=table.mean_log1p_difference.abs()).sort_values(
        ['q_BH_all_tested_contrasts','abs_effect','gene'],ascending=[True,False,True],kind='stable')
    top=ranked.groupby(['condition','label'],sort=False).head(10)
    top.to_csv(dest/'top10_descriptive_genes_per_contrast.csv',index=False)
    shown=ranked.gene.drop_duplicates().head(24).tolist()
    wide=table.pivot(index=['condition','label'],columns='gene',values='mean_log1p_difference').reindex(columns=shown)
    matrix=wide.to_numpy();limit=max(.1,float(np.nanmax(np.abs(matrix))))
    fig,ax=plt.subplots(figsize=(14,max(5,.32*len(wide))),layout='constrained')
    im=ax.imshow(np.ma.masked_invalid(matrix),aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit)
    ax.set_xticks(range(len(shown)),shown,rotation=65,ha='right',fontsize=8)
    ax.set_yticks(range(len(wide)),[('Fixed' if c=='fixed_glioma' else 'Selected')+' / '+lab for c,lab in wide.index],fontsize=8)
    ax.set_title('Unknown minus called within original cell class\nPatient-paired log-normalized expression; exploratory algorithm-status associations')
    fig.colorbar(im,ax=ax,label='Mean log1p-expression difference (not log2 fold change)')
    for ext in ['png','pdf']:fig.savefig(dest/f'Unknown_expression_contrasts.{ext}',dpi=210,bbox_inches='tight')
    plt.close(fig)
    note='''# Unknown expression contrasts: scope and assumptions

This exploratory analysis compares terminal Unknown with called cells within the
same author L1 cell class and sample. Both groups require at least20 cells. Sample
mean log1p(count/library total x10000) expression and detection fractions are
averaged within each patient; samples and cells are not independent replicates.
The matrix is the same upstream eligible-gene count matrix used by native R.
Normalization is checked against R's exported HVG expression for16 fixed-index
cells per sample. Missing sample-filtered genes are unavailable, not invented zeros.

The two conditions are fixed glioma markers and the training-patient-selected
marker on the same original HVG2000/UMAP/HDBSCAN partition. An eligible gene needs
at least6 paired patients and mean detection fraction>=0.1 in either status.
All class/condition support, including insufficient comparisons, is exported.
Patient means, detection fractions and every patient's effect are retained.

The primary descriptive test is an exact two-sided sign test on patient effects,
ignoring ties within1e-12. It tests sign balance among nonzero differences and
does not assume normality or symmetric effect magnitudes. Patient independence,
the frozen cohort and missing-gene availability delimit its interpretation.
Mean-effect percentile intervals use2000 patient resamples and retain observed
gene availability. The intervals estimate a different functional from the sign
test and are not used as its hypothesis test. BH correction is reported over all
tested genes/classes/conditions; within-contrast BH is secondary. The top-gene
figure is descriptive and does not imply that every displayed gene passes FDR.

Unknown is produced using the same expression data, so these associations are
partly induced by annotation and abstention rules. They are not independent
evidence of novel disease populations or causal mechanisms. Within-class pairing
reduces composition confounding but does not remove QC, depth or cell-state
differences. Inspect these results with the saved QC and marker-retention tables.
No doublet labels, transitional states or therapeutic targets are inferred here.
Recovery options to test later include broader supported marker/reference
coverage and calibrated abstention, with external validation before biological
interpretation; no new labels are assigned by this diagnostic.
'''
    (dest/'INTERPRETATION.md').write_text(note)
    write_json(dest/'manifest.json',dict(status='completed',exploratory=True,minimum_cells_per_status=MIN_CELLS,
        minimum_paired_patients=MIN_PATIENTS,minimum_mean_detection_fraction=.1,bootstrap_draws=BOOTSTRAPS,
        n_primary_samples=len(cohort),n_primary_patients=cohort.patient.nunique(),n_tested_gene_contrasts=len(table),
        test='Exact paired-patient sign test; global BH; patient-bootstrap intervals for mean effects',
        no_new_annotation_fits=True,not_independent_biological_discovery=True,
        source_manifests=sources,source_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest)

if __name__=='__main__':run()
