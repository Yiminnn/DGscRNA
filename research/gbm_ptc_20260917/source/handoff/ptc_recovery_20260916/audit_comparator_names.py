"""Map compact historical comparator spellings through the frozen native ontology."""
from pathlib import Path
import json
import os
import re
from ptc_common import BASE,RECOVERY,require_slurm,sha,utc,write_json
from label_rules import broad_lineage,simplify
from prepare_evaluation_reference import detection_metrics
require_slurm()
import pandas as pd
import numpy as np

def decode_name(name):
    if 'Î' in name or 'Ã' in name:
      try:return name.encode('latin1').decode('utf-8')
      except (UnicodeEncodeError,UnicodeDecodeError):pass
    return name

def token(name):
    name=re.sub(r'[^\w]+','',decode_name(name).lower(),flags=re.UNICODE).replace('_','')
    return name[:-1] if name.endswith('s') else name

ontology=pd.read_csv(BASE/'protocol/native_label_ontology.csv',keep_default_na=False)
lookup={}
for row in ontology.itertuples():
    key=token(row.general)
    if key in lookup:assert lookup[key]==row.broad,(key,lookup[key],row.broad)
    lookup[key]=row.broad

def mapped(name):
    general=simplify(decode_name(str(name)))
    key=token(general)
    if key in lookup:return lookup[key],'exact_frozen_ontology_name_after_format_normalization'
    return broad_lineage(general),'original_name_rule_fallback'

for name,want in [('Tcell','T'),('TregCells','T'),('Gammadelta(Î³Î´)Tcell','T'),
                  ('Naturalkillercell','NK'),('NaturalKillerT(NKT)cell','NKT'),
                  ('Plasmacytoiddendriticcell(pDC)','Myeloid')]:
    assert mapped(name)[0]==want,(name,mapped(name),want)
dest=BASE/'archived_comparators';dest.mkdir(exist_ok=True)
ref=pd.read_csv(BASE/'evaluation_reference/reference_cells.csv.gz',index_col=0,keep_default_na=False)
methods=['SCINA_archived','scCATCH_archived','scType_archived','SignacX_archived']
rows=[];dictionary=[]
for method in methods:
    native=ref[method+'_native']
    mapping={name:mapped(name) for name in native.unique()}
    labels=native.map({k:v[0] for k,v in mapping.items()})
    for name,(broad,rule) in mapping.items():
      dictionary.append(dict(method=method,native=name,decoded=decode_name(name),canonical_token=token(simplify(name)),
         broad=broad,rule=rule,prior_literal_rule=broad_lineage(name),n_cells=int(native.eq(name).sum())))
    for mapping_name,selected in [('strict_T',['T']),('legacy_broad_sensitivity',['T','NKT','NK','Lymphoid_ambiguous'])]:
      pred=labels.isin(selected)
      for tcr in ['TCR_cell_high_confidence_productive_TCR','TCR_any_filtered_contig','TCR_S3_supplied','TCR_paired_productive_TRA_TRB']:
        for sample,index in [('ALL',ref.index)]+[(s,g.index) for s,g in ref.groupby('sample')]:
          rows.append(dict(method=method,mapping=mapping_name,TCR_definition=tcr,sample=sample,
                 **detection_metrics(pred.loc[index],ref.loc[index,tcr])))
pd.DataFrame(dictionary).to_csv(dest/'format_normalized_name_ontology.csv',index=False)
pd.DataFrame(rows).to_csv(dest/'comparator_TCR_metrics.csv',index=False)
old=pd.read_csv(BASE/'evaluation_reference/historical_TCR_endpoint_replay.csv')
corrected=pd.concat([old[~old.method.isin(methods)],pd.DataFrame(rows)],ignore_index=True)
corrected.to_csv(dest/'historical_TCR_endpoint_replay_corrected_names.csv',index=False)
encoding=[]
for p in (RECOVERY/'s3_reconciliation').glob('*.SCINA_ct.mismatches.csv.gz'):
    d=pd.read_csv(p,index_col=0,keep_default_na=False)
    equal=d.S3.map(decode_name).eq(d.archive.map(decode_name))
    assert equal.all()
    encoding.append(dict(file=p.name,n_original_literal_mismatches=len(d),n_equal_after_UTF8_Latin1_repair=int(equal.sum())))
pd.DataFrame(encoding).to_csv(dest/'SCINA_623_encoding_parity.csv',index=False)
write_json(dest/'manifest.json',dict(status='completed',source_sha256=sha(Path(__file__)),
    ontology_sha256=sha(BASE/'protocol/native_label_ontology.csv'),job=os.environ['SLURM_JOB_ID'],completed_at=utc(),
    method='Normalize spacing,punctuation,terminal plural-s and recover evident UTF8/Latin1 mojibake; exact lookup through pre-frozen native ontology, otherwise original rule',
    mapping_does_not_use_TCR_outcomes=True,
    scope='Archived comparator spellings only; new DG-scRNA predictions and S2/S3 reference fields used by ongoing evaluation are unchanged',
    replaces='The initial historical comparator native-name metrics had compact Tcell/Thelper names incorrectly outside the T rule; they are diagnostic drafts and must not be used to rank methods.',
    source_flags='Literal S3 binary flags retained as a separate original-paper sensitivity'))
(dest/'COMPLETE').write_text(sha(dest/'manifest.json')+'\n')
print('Archived compact-name ontology corrected; all623 SCINA differences are encoding-only')
