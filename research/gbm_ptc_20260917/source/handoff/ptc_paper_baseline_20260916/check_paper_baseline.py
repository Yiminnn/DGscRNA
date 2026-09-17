"""Recompute the literal paper/S3 endpoint and trace group-specific output provenance.

This is a baseline reconciliation job, not an ablation or a search for a metric
definition that agrees with the manuscript. All source flags remain unchanged.
"""
from pathlib import Path
import hashlib
import json
import os
import re
from zipfile import ZipFile
import xml.etree.ElementTree as ET

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
BASE = ROOT / 'results/hvg_ptc_20260916_v1'
REC = BASE / 'ptc_recovery'
OUT = BASE / 'ptc_paper_baseline'


def run():
    assert os.environ.get('SLURM_JOB_ID'), 'Use SLURM'
    import numpy as np
    import pandas as pd
    from sklearn.metrics import f1_score, accuracy_score, roc_auc_score, confusion_matrix
    OUT.mkdir(exist_ok=True)
    archive = REC/'archive/tcr'
    metadata = pd.read_csv(archive/'rawdata/metadata.txt', sep='\t')
    samplemap = dict(zip(metadata.sc_ID.astype(str), metadata.Sample))
    sources = [archive/'rawdata/data_with_validation.csv',
               archive/'rawdata/data_with_validation+3_cell_types 2.csv',
               ROOT/'paper/submission_v16/Supplementary table S3_Method comparison in cell type annotation.xlsx',
               REC/'inventory_workspace/object_01.meta.csv',
               archive/'ptc_val/scripts/DGscRNA-Share/R/source.R',
               archive/'ptc_val/scripts/DGscRNA-Vignette.rmd']

    def canonical(frame, idfield):
        ids = frame[idfield].astype(str).str.extract(r'([ACGT]{12,})', expand=False)
        frame = frame.copy()
        frame['cell_id'] = frame['sample.name'].map(samplemap) + '_' + ids
        assert frame.cell_id.notna().all() and frame.cell_id.is_unique
        frame['sample'] = frame['sample.name'].map(samplemap)
        frame['analysis_group'] = np.where(frame['sample'].str.startswith(('MT-', 'N-')), 'NMT', 'TTU')
        frame['paper_tissue'] = np.select([frame['sample'].str.startswith('N-'),frame['sample'].str.startswith('MT-')],
                                         ['Normal adjacent','Metastasis'], default='Primary tumor')
        return frame.set_index('cell_id')

    old = pd.read_csv(sources[0], keep_default_na=False)
    old = canonical(old, old.columns[0])
    extended = pd.read_csv(sources[1], keep_default_na=False)
    extended = canonical(extended, 'X')
    s3 = pd.read_excel(sources[2], header=5, keep_default_na=False)
    s3.columns = s3.columns.astype(str).str.strip()
    s3 = s3[s3['sample.name'].isin(samplemap)]
    s3 = canonical(s3, 'Sample_ID')
    assert len(old)==len(s3)==len(extended)==92404
    assert set(old.index)==set(s3.index)==set(extended.index)
    s3=s3.loc[old.index]; extended=extended.loc[old.index]
    parity=[]
    for target,source in [('DGscRNA_prediction_cell_type','Final_DGCyTOF'),
                          ('DGscRNA_prediction_detailed_T_cell','Final_DGCyTOF_General_w_T_cell'),
                          ('DGscRNA_ct_T_cells','T_cells'),('T_cell','validation_t_cell')]:
        a=s3[target].astype(str);b=old[source].astype(str)
        if target in ['DGscRNA_ct_T_cells','T_cell']:
            a=pd.to_numeric(s3[target]);b=pd.to_numeric(old[source])
        parity.append(dict(S3_column=target,archive_column=source,n=len(a),n_equal=int(a.eq(b).sum()),
                           n_mismatch=int(a.ne(b).sum())))
    pd.DataFrame(parity).to_csv(OUT/'S3_literal_output_parity.csv',index=False)

    # Source-defined broad historical T labels; do not modify to optimize F1.
    vignette=sources[-1].read_text()
    section=vignette.split('T_cell_types = c(',1)[1].split('\n\nmetadata',1)[0]
    historical_T=re.findall(r"'([^']*)'",section)
    (OUT/'historical_T_names_from_vignette.json').write_text(json.dumps(historical_T,indent=2)+'\n')
    predictions = {'DG_S3_literal_flag':pd.to_numeric(s3.DGscRNA_ct_T_cells).astype(int),
        'DG_archived_general_with_vignette_T_list':old.Final_DGCyTOF_General.isin(historical_T).astype(int),
        'scCATCH_S3_literal_flag':pd.to_numeric(s3.scCATCH_ct_T_cells).astype(int),
        'SCINA_S3_literal_flag':pd.to_numeric(s3.SCINA_ct_T_cells).astype(int),
        'scType_S3_literal_flag':pd.to_numeric(s3.scType_ct_T_cells).astype(int),
        'SignacX_S3_literal_flag':pd.to_numeric(s3.signacx_ct_T_cells).astype(int)}
    truth=pd.to_numeric(s3.T_cell).astype(int)
    assert set(truth)=={0,1}
    scopes=[('Overall',old.index)]+list(old.groupby('paper_tissue').groups.items())+list(old.groupby('analysis_group').groups.items())+list(old.groupby('sample').groups.items())
    rows=[]
    for name,pred in predictions.items():
        assert set(pred)<= {0,1}
        for scope,ids in scopes:
            y=truth.loc[ids];p=pred.loc[ids]
            tn,fp,fn,tp=confusion_matrix(y,p,labels=[0,1]).ravel()
            rows.append(dict(method=name,scope=scope,n=len(ids),TCR_positive=int(y.sum()),
                predicted_positive=int(p.sum()),TP=int(tp),FP=int(fp),FN=int(fn),TN=int(tn),
                F1_T_positive=f1_score(y,p,pos_label=1,zero_division=0),
                F1_nonT_positive_diagnostic=f1_score(y,p,pos_label=0,zero_division=0),
                F1_macro_diagnostic=f1_score(y,p,average='macro',zero_division=0),
                F1_weighted_diagnostic=f1_score(y,p,average='weighted',zero_division=0),
                accuracy=accuracy_score(y,p),AUC_from_binary_calls=roc_auc_score(y,p)))
    metrics=pd.DataFrame(rows)
    metrics.to_csv(OUT/'paper_S3_literal_metrics.csv',index=False)

    # Read targets directly from the submitted document, not from an earlier memo.
    document=ROOT/'paper/submission_v16/DG_scRNA_04232026_V16_cell_report.docx'
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    with ZipFile(document) as z:xml=ET.fromstring(z.read('word/document.xml'))
    tables=[]
    for t in xml.findall('.//w:tbl',ns):
        rows=[]
        for r in t.findall('./w:tr',ns):
            rows.append([''.join(x.text or '' for x in c.findall('.//w:t',ns)) for c in r.findall('./w:tc',ns)])
        tables.append(rows)
    table=next(t for t in tables if any('Datasets' in r and 'DG-scRNA' in r for r in t))
    (OUT/'paper_Table2_literal_cells.json').write_text(json.dumps(table,indent=2)+'\n')
    targets=[];scope=None
    method_columns={'scCATCH':2,'SCINA':3,'scType':4,'DG':6}
    method_names={'scCATCH':'scCATCH_S3_literal_flag','SCINA':'SCINA_S3_literal_flag','scType':'scType_S3_literal_flag','DG':'DG_S3_literal_flag'}
    metric_columns={'F1 score':'F1_T_positive','Accuracy':'accuracy','AUC-ROC':'AUC_from_binary_calls'}
    for row in table:
        if len(row)!=7 or row[1] not in metric_columns: continue
        if row[0]:scope=row[0]
        for method,col in method_columns.items():
            target=float(row[col]);key=metric_columns[row[1]]
            value=float(metrics[(metrics.method==method_names[method]) & (metrics.scope==scope)][key].iloc[0])
            targets.append(dict(method=method,scope=scope,metric=row[1],paper=target,
                literal_S3_recalculated=value,difference=value-target,
                agrees_at_reported_4_decimals=round(value,4)==target,
                calculation='T-cell positive F1/accuracy; AUC uses supplied binary predictions, historical probabilities unavailable'))
    pd.DataFrame(targets).to_csv(OUT/'paper_Table2_target_comparison.csv',index=False)

    # Check branch provenance using only saved outputs, never fit to reference labels.
    branchcols=[c for c in extended if c.startswith('hdbscan.UMAP_clusters_')]
    branchchecks=[]
    for group,ids in old.groupby('analysis_group').groups.items():
        for c in branchcols:
            equal=old.loc[ids,'Final_DGCyTOF_General'].eq(extended.loc[ids,c])
            branchchecks.append(dict(group=group,archived_branch=c,n=len(ids),n_equal=int(equal.sum()),
                                     n_mismatch=int((~equal).sum())))
    pd.DataFrame(branchchecks).to_csv(OUT/'archived_selected_branch_checks.csv',index=False)
    inventories=[]
    for field in ['Final_DGCyTOF','Final_DGCyTOF_General','Final_DGCyTOF_General_w_T_cell','pred_cell_types']:
        for (group,label),n in old.groupby(['analysis_group',field]).size().items():
            inventories.append(dict(group=group,field=field,label=label,n_cells=int(n)))
    pd.DataFrame(inventories).to_csv(OUT/'original_final_labels_by_group.csv',index=False)

    # The saved S2 object need not be the later paper's final selection.
    saved=pd.read_csv(sources[3],keep_default_na=False)
    saved=canonical(saved,saved.columns[0]).loc[old.index]
    cluster_match=old.seurat_clusters.astype(str).eq(saved.seurat_clusters.astype(str))
    pd.DataFrame({'cell_id':old.index,'paper_seurat':old.seurat_clusters.values,
                  'checkpoint_seurat':saved.seurat_clusters.values,'exact':cluster_match.values}).to_csv(
                      OUT/'checkpoint_paper_seurat_partition_parity.csv.gz',index=False)
    pd.DataFrame(dict(cell_id=old.index,sample=old['sample'],group=old.analysis_group,
                     seurat_clusters=old.seurat_clusters,hdbscan_UMAP_clusters=old['hdbscan.UMAP_clusters'],
                     paper_final_native=old.Final_DGCyTOF,paper_final_general=old.Final_DGCyTOF_General,
                     paper_T_flag=truth,DG_T_flag=predictions['DG_S3_literal_flag'])).to_csv(
                         OUT/'paper_baseline_reference.csv.gz',index=False)
    manifest=dict(status='baseline_gate_under_reconciliation',job=os.environ['SLURM_JOB_ID'],
        n_cells=len(old),checkpoint_seurat_exact_cells=int(cluster_match.sum()),
        paper_Table2_rows_recomputed=len(targets),
        paper_Table2_rows_matching_4_decimals=sum(x['agrees_at_reported_4_decimals'] for x in targets),
        source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources+[document,Path(__file__)]})
    (OUT/'literal_reconciliation_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(metrics[(metrics.method=='DG_S3_literal_flag') & metrics.scope.isin(['Overall','Normal adjacent','Primary tumor','Metastasis'])].to_string(index=False))
    print(pd.DataFrame(targets).query("method=='DG'").to_string(index=False))
    print('Saved Seurat partition exact:',int(cluster_match.sum()),'/',len(old))
    print(pd.DataFrame(branchchecks).to_string(index=False))


if __name__=='__main__':
    run()
