"""Evaluation-only joins and historical endpoint replay, after ontology is frozen."""
from pathlib import Path
import json
import os
from ptc_common import ROOT,BASE,RECOVERY,require_slurm,sha,write_json,utc,GROUPS
from label_rules import simplify,broad_lineage,strict_T,legacy_broad_T_sensitivity

def detection_metrics(pred,detected):
    import numpy as np
    pred=np.asarray(pred,dtype=bool);detected=np.asarray(detected,dtype=bool)
    tp=int((pred&detected).sum());fp=int((pred&~detected).sum())
    fn=int((~pred&detected).sum());tn=int((~pred&~detected).sum())
    return dict(n=len(pred),n_predicted_T=int(pred.sum()),n_detected_TCR=int(detected.sum()),
      detected_and_predicted=tp,detected_not_predicted=fn,predicted_no_detection=fp,neither=tn,
      TCR_positive_recall=tp/(tp+fn) if tp+fn else None,
      TCR_detection_yield_in_predicted_T=tp/(tp+fp) if tp+fp else None,
      apparent_TCR_binary_F1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0.0)

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    dest=BASE/'evaluation_reference';dest.mkdir(exist_ok=True)
    frozen=json.loads((BASE/'protocol/label_rules_manifest.json').read_text())
    assert frozen['rule_source_sha256']==sha(Path(__file__).with_name('label_rules.py'))
    paths=dict(QC=RECOVERY/'fixed_QC_inputs/cells_QC_sample_patient.csv',
      S2=RECOVERY/'table_reconciliation/S2_as_supplied.csv.gz',
      S3=RECOVERY/'s3_reconciliation/S3_as_supplied_canonical_ids.csv.gz',
      TCR=RECOVERY/'tcr_audit/all_S2_TCR_detection_rules.csv.gz')
    qc=pd.read_csv(paths['QC'],index_col=0)
    s2=pd.read_csv(paths['S2'],index_col=0,keep_default_na=False)
    s3=pd.read_csv(paths['S3'],index_col=0,keep_default_na=False)
    tcr=pd.read_csv(paths['TCR'],index_col=0)
    for d in [qc,s2,s3,tcr]:assert d.index.is_unique and set(d.index)==set(qc.index)
    s2=s2.loc[qc.index];s3=s3.loc[qc.index];tcr=tcr.loc[qc.index]
    ref=qc.copy();ref.index.name='cell_id'
    ref['sample']=ref.sample_id.astype(str);ref['patient']=ref.patient_id.astype(str)
    ref['group']=ref['sample'].map({s:g for g,ss in GROUPS.items() for s in ss})
    ref['tissue_name_code']=ref['sample'].str.split('-').str[0]
    for col in tcr:
      if col not in ['sample','patient']:ref['TCR_'+col]=tcr[col]
    ref['TCR_S3_supplied']=pd.to_numeric(s3.T_cell,errors='raise').astype(bool)
    annotations={'S2_initial':s2.Cell_type_annotation_CellMarker2_0 if 'Cell_type_annotation_CellMarker2_0' in s2 else s2['Cell_type_annotation_CellMarker2.0'],
       'S2_terminal':s2['DG_scRNA_Finalized_Cell_Types'],
       'S3_terminal':s3['DGscRNA_prediction_cell_type'],
       'S3_detailed_T':s3['DGscRNA_prediction_detailed_T_cell'],
       'SCINA_archived':s3['SCINA_ct'],'scCATCH_archived':s3['scCATCH'],
       'scType_archived':s3['scType'],'SignacX_archived':s3['signacx']}
    maprows=[]
    for method,values in annotations.items():
      ref[method+'_native']=values
      mapping={name:broad_lineage(name) for name in values.unique()}
      ref[method+'_broad']=values.map(mapping)
      for name in values.unique():maprows.append(dict(method=method,native=name,general=simplify(name),
                    broad=mapping[name],strict_T=strict_T(name),legacy_broad_T=legacy_broad_T_sensitivity(name)))
    source_flags={'DGscRNA_S3_supplied_T_flag':'DGscRNA_ct_T_cells',
      'SCINA_S3_supplied_T_flag':'SCINA_ct_T_cells','scCATCH_S3_supplied_T_flag':'scCATCH_ct_T_cells',
      'scType_S3_supplied_T_flag':'scType_ct_T_cells','SignacX_S3_supplied_T_flag':'signacx_ct_T_cells'}
    for name,column in source_flags.items():ref[name]=pd.to_numeric(s3[column],errors='raise').astype(bool)
    modules=[]
    for group in GROUPS:
      assert (BASE/f'biology/{group}.COMPLETE').exists()
      modules.append(pd.read_csv(BASE/f'biology/{group}.module_scores.csv.gz').set_index('cell_id'))
    module=pd.concat(modules).loc[ref.index]
    assert module.index.is_unique and not module.isna().any().any()
    ref=ref.join(module)
    assert len(ref)==92404 and ref['group'].notna().all() and ref.patient.nunique()==4
    ref.to_csv(dest/'reference_cells.csv.gz')
    pd.DataFrame(maprows).to_csv(dest/'historical_name_ontology.csv',index=False)
    records=[]
    for method,values in annotations.items():
      for mapping,rule in [('strict_T',strict_T),('legacy_broad_sensitivity',legacy_broad_T_sensitivity)]:
        pred=values.map({v:rule(v) for v in values.unique()})
        for rule_name in ['TCR_cell_high_confidence_productive_TCR','TCR_any_filtered_contig','TCR_S3_supplied',
                          'TCR_paired_productive_TRA_TRB']:
          for sample,index in [('ALL',ref.index)]+[(s,g.index) for s,g in ref.groupby('sample')]:
            records.append(dict(method=method,mapping=mapping,TCR_definition=rule_name,sample=sample,
              **detection_metrics(pred.loc[index],ref.loc[index,rule_name])))
    for method in source_flags:
      for rule_name in ['TCR_cell_high_confidence_productive_TCR','TCR_any_filtered_contig','TCR_S3_supplied']:
        for sample,index in [('ALL',ref.index)]+[(s,g.index) for s,g in ref.groupby('sample')]:
          records.append(dict(method=method,mapping='literal_source_binary_flag',TCR_definition=rule_name,sample=sample,
               **detection_metrics(ref.loc[index,method],ref.loc[index,rule_name])))
    pd.DataFrame(records).to_csv(dest/'historical_TCR_endpoint_replay.csv',index=False)
    libs=json.loads((RECOVERY/'inventory_r/object_01_full_marker_symbol.content.json').read_text())
    retention=[]
    for group,samples in GROUPS.items():
      scoring=set((BASE/f'prepared/{group}/RNA_genes.txt').read_text().splitlines())
      for sample in samples:
        eligible=(BASE/f'samples/{sample}/geometry_genes.txt').read_text().splitlines()
        ranked=(BASE/f'samples/{sample}/HVG_rank_top5000.txt').read_text().splitlines()
        for feature in ['all','hvg500','hvg1000','hvg2000','hvg3000','hvg5000']:
          selected=set(eligible if feature=='all' else ranked[:int(feature[3:])])
          for library,panels in libs.items():
            for panel,genes in panels.items():
              genes=[genes] if isinstance(genes,str) else genes
              original=set(genes)
              retention.append(dict(sample=sample,group=group,feature=feature,library=library,panel=panel,
                 full_panel_denominator=len(genes),unique_marker_symbols=len(original),
                 scoring_retained=len(original&scoring),geometry_retained=len(original&selected),
                 scoring_retained_but_geometry_lost=len((original&scoring)-selected),n_geometry_features=len(selected)))
    pd.DataFrame(retention).to_csv(dest/'geometry_marker_retention.csv.gz',index=False)
    pd.DataFrame([dict(sample=s,patient=g.patient.iloc[0],group=g.group.iloc[0],n_cells=len(g),
       TCR_productive_positive=int(g.TCR_cell_high_confidence_productive_TCR.sum()),
       any_contig_positive=int(g.TCR_any_filtered_contig.sum()),S3_positive=int(g.TCR_S3_supplied.sum()),
       median_nCount_RNA=float(g.nCount_RNA.median()),median_nFeature_RNA=float(g.nFeature_RNA.median()),
       median_percent_mt=float(g['percent.mt'].median())) for s,g in ref.groupby('sample')]).to_csv(dest/'cohort_summary.csv',index=False)
    write_json(dest/'manifest.json',dict(status='completed',source_sha256=sha(Path(__file__)),
       inputs={k:sha(v) for k,v in paths.items()},label_rules_sha256=sha(BASE/'protocol/label_rules_manifest.json'),
       n_cells=len(ref),n_patients=ref.patient.nunique(),completed_at=utc(),job=os.environ['SLURM_JOB_ID'],
       primary_TCR_rule='cell_high_confidence_productive_TCR',
       binary_F1_interpretation='Agreement with detection indicator; assumes undetected=negative only for this descriptive assay metric, not biological accuracy',
       competitor_scope='Recovered original outputs; original pooling/reference/input conditions differ from new controlled branches'))
    (dest/'COMPLETE').write_text(sha(dest/'manifest.json')+'\n')
    print('Evaluation reference, historical endpoints and marker retention complete')

if __name__=='__main__':run()
