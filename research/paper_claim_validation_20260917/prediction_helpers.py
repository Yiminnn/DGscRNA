"""Read frozen native predictions and map terms without consulting evaluation labels."""
import json
from functools import lru_cache
from common import OUT, L1, require_slurm, sha

UNKNOWN = {'unknown','Unknown','Undecided','','NA','Noise'}

def map_labels(method, library, labels):
    require_slurm()
    import numpy as np
    import pandas as pd
    if method=='SingleR':
        lookup={v:v for v in L1}
    elif method=='scDeepSort':
        from build_markers import semantic
        lookup={v:semantic(v) for v in set(labels)}
    else:
        table=pd.read_csv(OUT/'markers/panel_L1_mapping.csv',dtype=str,keep_default_na=False)
        table=table[table.library==library];assert len(table)
        lookup=dict(zip(table.panel,table.L1))
    @lru_cache(None)
    def split(value):
        if value in lookup:return [(value,)]
        found=[]
        for panel in lookup:
            if value.startswith(panel+', '):found.extend((panel,)+rest for rest in split(value[len(panel)+2:]))
        return found
    def one(value):
        if value in lookup:return lookup[value]
        if value in UNKNOWN:return 'Unknown'
        if method=='scCATCH':
            possibilities=split(value)
            parents={lookup[p] for row in possibilities for p in row}
            if possibilities and len(parents)==1:return next(iter(parents))
        return 'UNMAPPABLE'
    mapping={v:one(v) for v in set(labels)}
    return np.asarray([mapping[v] for v in labels])

def read_native(method,sample,budget,route,library,cutoff):
    require_slurm()
    import pandas as pd
    if method=='DG-scRNA':
        path=OUT/'GBM'/sample/budget/route
        manifest=json.loads((path/'score_manifest.json').read_text())
        matches=[a for a,p in manifest['arms'].items() if p['library']==library and str(p['cutoff'])==str(cutoff)]
        assert len(matches)==1,(method,sample,library,cutoff)
        source=path/'terminal'/matches[0]/'predictions.csv.gz'
        data=pd.read_csv(source,dtype=str,keep_default_na=False)
        return data.rename(columns={'final090':'prediction'}),source
    path=OUT/'comparators'/method/sample
    if method=='scType':
        manifest=json.loads((path/'cohort_manifest.json').read_text());source=path/'cohort_predictions.csv.gz'
    elif method=='scCATCH':
        path=path/'hvg2000/UMAP2_HDBSCAN_R'
        manifest=json.loads((path/'cohort_manifest.json').read_text());source=path/'cohort_predictions.csv.gz'
    elif method=='SCINA':
        libraries=list(json.loads((OUT/'markers/libraries.json').read_text()))
        path=path/f'L{libraries.index(library):02d}'
        manifest=json.loads((path/'manifest.json').read_text());source=path/'predictions.csv.gz'
    elif method in ['SingleR','scDeepSort']:
        source=path/'predictions.csv.gz'
        data=pd.read_csv(source,dtype=str,keep_default_na=False)
        return data[['cell_id',cutoff]].rename(columns={cutoff:'prediction'}),source
    else:raise ValueError(method)
    matches=[a for a,p in manifest['arms'].items() if p['library']==library and str(p['cutoff'])==str(cutoff)
             and p['budget']==budget and p['route']==route]
    assert len(matches)==1,(method,sample,budget,route,library,cutoff,matches)
    assert sha(source)==manifest['predictions_sha256']
    data=pd.read_csv(source,dtype=str,keep_default_na=False)
    return data[['cell_id',matches[0]]].rename(columns={matches[0]:'prediction'}),source
