"""All-cell comparator scoring with the same frozen native-panel ontology mapping."""
import json
import os
import sys
from pathlib import Path
from common import OUT, L1, require_slurm, checked, complete, sha, write_json, utc

def run(method,sample):
    require_slurm()
    import numpy as np
    import pandas as pd
    from sklearn.metrics import precision_recall_fscore_support
    root=OUT/'comparators'/method/sample;dest=root/'evaluation'
    if checked(dest):return
    dest.mkdir(parents=True,exist_ok=True)
    inputs=json.loads((OUT/'inputs'/sample/'input_manifest.json').read_text())
    tp=OUT/'evaluation_inputs'/sample/'truth.csv.gz'
    assert sha(tp)==inputs['evaluation_files']['truth.csv.gz']
    truth=pd.read_csv(tp,dtype=str,keep_default_na=False);y=truth.L1.to_numpy()
    mapping=pd.read_csv(OUT/'markers/panel_L1_mapping.csv',dtype=str,keep_default_na=False)
    maps={lib:dict(zip(g.panel,g.L1)) for lib,g in mapping.groupby('library')}
    if method=='scType':
        sources=[(root,'cohort_manifest.json','COHORT_COMPLETE','cohort_predictions.csv.gz')]
    elif method=='SCINA':
        sources=[(root/f'L{i:02d}','manifest.json','COMPLETE','predictions.csv.gz') for i in range(16)]
    elif method in ['SingleR','scDeepSort']:
        sources=[(root,'manifest.json','COMPLETE','predictions.csv.gz')]
    elif method=='scCATCH':
        sources=[(root/'hvg2000/UMAP2_HDBSCAN_R','cohort_manifest.json','COHORT_COMPLETE','cohort_predictions.csv.gz')]
    else:raise ValueError(method)
    rows=[];classes=[];confusions=[];files={}
    for src,manifest,flag,predfile in sources:
        assert checked(src,manifest,flag),str(src)
        m=json.loads((src/manifest).read_text());assert m['status']=='completed'
        if method=='SingleR':
            m['arms']={key:dict(library='labelled_training_patients',budget='shared_reference_genes',route='cell_level',cutoff=key)
                       for key in ['default','unpruned']}
            maps['labelled_training_patients']={v:v for v in L1}
        assert sha(src/predfile)==m['predictions_sha256']
        pred=pd.read_csv(src/predfile,dtype=str,keep_default_na=False)
        if method=='scDeepSort':
            from build_markers import semantic
            m['arms']={'default':dict(library='published_human_Brain_atlas',budget='pretrained_genes',route='GNN',cutoff='default')}
            maps['published_human_Brain_atlas']={v:semantic(v) for v in set(pred.default)}
            maps['published_human_Brain_atlas']['Unknown']='Unknown'
        assert np.array_equal(pred.cell_id,truth.cell_id)
        files[str(src/manifest)]=sha(src/manifest)
        for aid,arm in m['arms'].items():
            assert arm.get('status','completed')!='implementation_error'
            native=pred[aid].to_numpy()
            abstain=np.isin(native,['unknown','Unknown','Undecided','','NA'])
            lookup=maps[arm['library']]
            def mapped(value,is_unknown):
                if value in lookup:return lookup[value]
                if is_unknown:return 'Unknown'
                if method=='scCATCH':
                    # The official tool returns all tied native panels. A tie may
                    # collapse to one L1 parent only if every exact panel agrees.
                    from functools import lru_cache
                    @lru_cache(None)
                    def split(s):
                        if s in lookup:return [[s]]
                        found=[]
                        for panel in lookup:
                            if s.startswith(panel+', '):found.extend([[panel]+r for r in split(s[len(panel)+2:])])
                        return found
                    possibilities=split(value)
                    if possibilities:
                        parents={lookup[v] for row in possibilities for v in row}
                        if len(parents)==1:return next(iter(parents))
                return 'UNMAPPABLE'
            value_map={v:mapped(v,v in ['unknown','Unknown','Undecided','','NA']) for v in set(native)}
            p=np.asarray([value_map[v] for v in native])
            pr,re,f1,su=precision_recall_fscore_support(y,p,labels=L1,zero_division=0)
            ctx=dict(method=method,sample=sample,patient=inputs['patient'],primary=inputs['primary'],
                budget=arm['budget'],route=arm['route'],library=arm['library'],cutoff=arm['cutoff'],
                arm_id=aid,status=arm.get('status','completed'),n_cells=len(y),
                reference_overlap=arm['library'] in ['CARE_TME','BrainAtlas112','UNION_all'])
            rows.append(dict(**ctx,macroF1_present=float(f1[su>0].mean()),macroF1_fixed11=float(f1.mean()),
                accuracy=float((p==y).mean()),coverage=float((~abstain).mean()),unknown_rate=float(abstain.mean()),
                mapped_coverage=float(np.isin(p,L1).mean()),off_vocabulary_rate=float((~np.isin(p,L1)&~abstain).mean())))
            for i,label in enumerate(L1):classes.append(dict(**ctx,label=label,precision=float(pr[i]),recall=float(re[i]),F1=float(f1[i]),support=int(su[i])))
            for r in pd.DataFrame({'truth':y,'prediction':p}).value_counts().reset_index(name='n').to_dict('records'):
                confusions.append(dict(**ctx,**r))
    pd.DataFrame(rows).to_csv(dest/'metrics.csv',index=False)
    pd.DataFrame(classes).to_csv(dest/'per_class.csv.gz',index=False,compression='gzip')
    pd.DataFrame(confusions).to_csv(dest/'confusions.csv.gz',index=False,compression='gzip')
    write_json(dest/'manifest.json',dict(status='completed',method=method,sample=sample,n_conditions=len(rows),
        full_cell_denominator=True,unknown_and_unmapped_count_as_errors=True,input_manifests=files,
        files={p.name:sha(p) for p in dest.iterdir() if p.suffix in ['.csv','.gz']},truth_sha256=sha(tp),
        job=os.environ['SLURM_JOB_ID'],source_sha256=sha(__file__),completed_at=utc()))
    complete(dest)

if __name__=='__main__':run(*sys.argv[1:])
