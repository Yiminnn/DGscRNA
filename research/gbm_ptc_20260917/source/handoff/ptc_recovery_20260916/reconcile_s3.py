#!/usr/bin/env python3
"""Reconcile S3 to archived terminal exports and preserve all S2/S3 differences."""
from pathlib import Path
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'hvg_ptc_20260916'))
from common import ROOT,OUT,require_slurm,sha,utc,write_json,runtime_record


def run():
    require_slurm()
    import pandas as pd
    base=OUT/'ptc_recovery';archive=base/'archive/tcr'
    assert (base/'tcr_audit/COMPLETE').exists()
    dest=base/'s3_reconciliation';dest.mkdir(exist_ok=True);start=time.monotonic()
    meta=pd.read_csv(archive/'rawdata/metadata.txt',sep='\t').set_index('Sample')
    samplemap=dict(zip(meta.sc_ID.astype(str),meta.index))
    book=ROOT/'paper/submission_v16/Supplementary table S3_Method comparison in cell type annotation.xlsx'
    s3=pd.read_excel(book,header=5)
    s3.columns=s3.columns.astype(str).str.strip()
    s3=s3.loc[s3['sample.name'].notna()].copy()
    s3['canonical_cell_id']=s3['sample.name'].map(samplemap)+'_'+s3.Sample_ID.str.extract(r'([ACGT]{12,})',expand=False)
    assert s3.canonical_cell_id.notna().all() and s3.canonical_cell_id.is_unique
    s3=s3.set_index('canonical_cell_id')
    s3.to_csv(dest/'S3_as_supplied_canonical_ids.csv.gz')
    rules=pd.read_csv(base/'tcr_audit/all_S2_TCR_detection_rules.csv.gz',index_col=0)
    assert rules.index.is_unique
    common=rules.index.intersection(s3.index)
    diff=rules.loc[common].copy()
    diff['S3_T_cell']=pd.to_numeric(s3.loc[common,'T_cell'],errors='raise').astype(int)
    diff['S2_S3_agree']=diff.S2_is_T_cell_Real.eq(diff.S3_T_cell)
    diff['S3_DGscRNA_prediction']=s3.loc[common,'DGscRNA_prediction_cell_type']
    diff['S3_DGscRNA_detailed_T_cell']=s3.loc[common,'DGscRNA_prediction_detailed_T_cell']
    diff.loc[~diff.S2_S3_agree].to_csv(dest/'all_S2_S3_TCR_disagreements.csv.gz')
    # Original exports are joined by sample + barcode, with no relabelling.
    correspondence={
        'DGscRNA_prediction_cell_type':'Final_DGCyTOF',
        'DGscRNA_prediction_detailed_T_cell':'Final_DGCyTOF_General_w_T_cell',
        'T_cell':'validation_t_cell',
        'SCINA_ct':'SCINA_ct_General','SCINA_ct_T_cells':'SCINA_ct_T_cells',
        'scCATCH':'scCATCH_ct_General','scCATCH_ct_T_cells':'scCATCH_ct_T_cells',
        'scType':'scType_ct_General','scType_ct_T_cells':'scType_ct_T_cells',
        'signacx':'signacx_ct','signacx_ct_T_cells':'signacx_ct_T_cells',
        'DGscRNA':'pred_cell_types','DGscRNA_ct_T_cells':'T_cells',
        'normal_or_cancer_General':'normal_or_cancer_General',
        'tissue_sites_ls':'tissue_sites_ls'}
    comparisons=[];provenance=[]
    for name in ['data_with_validation.csv','data_with_validation+3_cell_types 2.csv']:
        path=archive/'rawdata'/name
        d=pd.read_csv(path,index_col=0)
        d.index=d['sample.name'].map(samplemap)+'_'+d.index.to_series().str.extract(r'([ACGT]{12,})',expand=False)
        assert d.index.notna().all() and d.index.is_unique
        joined=s3.index.intersection(d.index)
        provenance.append(dict(source=str(path),sha256=sha(path),n_rows=len(d),n_joined=len(joined)))
        for target,source in correspondence.items():
            a=s3.loc[joined,target];b=d.loc[joined,source]
            numeric=target=='T_cell' or target.endswith('_T_cells')
            if numeric:a=pd.to_numeric(a);b=pd.to_numeric(b)
            else:a=a.astype('string');b=b.astype('string')
            equal=(a==b).fillna(False)
            comparisons.append(dict(source_file=name,S3_column=target,archive_column=source,
                n_joined=len(joined),n_equal=int(equal.sum()),n_mismatch=int((~equal).sum()),
                all_S3_exact=bool(equal.all() and len(joined)==len(s3))))
            if not equal.all():
                pd.DataFrame({'S3':a[~equal],'archive':b[~equal]}).to_csv(
                    dest/f'{path.stem}.{target}.mismatches.csv.gz')
    pd.DataFrame(comparisons).to_csv(dest/'S3_to_archive_column_agreement.csv',index=False)
    summary=[]
    for sample,g in diff.groupby('sample',sort=False):
        summary.append(dict(sample=sample,patient=g.patient.iloc[0],n_joined=len(g),
            S2_TCR_positive=int(g.S2_is_T_cell_Real.sum()),S3_TCR_positive=int(g.S3_T_cell.sum()),
            S2_positive_S3_negative=int((g.S2_is_T_cell_Real.eq(1)&g.S3_T_cell.eq(0)).sum()),
            S2_negative_S3_positive=int((g.S2_is_T_cell_Real.eq(0)&g.S3_T_cell.eq(1)).sum())))
    pd.DataFrame(summary).to_csv(dest/'S2_S3_detection_difference_by_sample.csv',index=False)
    # Barcode-only matching is diagnosed separately, never used as a repair.
    all_tcr_barcodes=set()
    for p in (base/'tcr_audit').glob('*.original_TCR_barcode_rules.csv.gz'):
        ids=pd.read_csv(p,usecols=['cell_id']).cell_id.str.split('_').str[-1]
        all_tcr_barcodes.update(ids)
    diff['barcode_seen_in_any_sample_TCR']=diff.index.to_series().str.split('_').str[-1].isin(all_tcr_barcodes)
    diff.to_csv(dest/'S2_S3_TCR_full_join.csv.gz')
    mismatch=diff.loc[~diff.S2_S3_agree]
    write_json(dest/'manifest.json',dict(status='completed_discrepancies_retained',timestamp=utc(),
        n_S2=len(rules),n_S3=len(s3),n_joined=len(common),n_S2_S3_disagreements=len(mismatch),
        mismatch_barcodes_seen_in_another_sample_TCR=int(mismatch.barcode_seen_in_any_sample_TCR.sum()),
        S3_sha256=sha(book),source_exports=provenance,source_sha256=sha(Path(__file__)),
        seconds=time.monotonic()-start,**runtime_record()))
    (dest/'COMPLETE').write_text(sha(dest/'manifest.json')+'\n')


if __name__=='__main__':run()
