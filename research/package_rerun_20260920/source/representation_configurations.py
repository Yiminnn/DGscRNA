"""Verbatim finite scientific configuration function from the archived controller."""
LIBRARIES=['CM2_glioma_other','CM2_primary_all_context']

def configurations(sample,budget,include_defaults=False):
    result=[]
    def add(space,method,minpts=50,res=.5,seed=42,kind='parameter'):
        name=f'{space}_{method}_minPts{minpts}_r{res}_seed{seed}'
        result.append(dict(name=name,sample=sample,budget=budget,space=space,method=method,
            minPts=minpts,resolution=res,embedding_seed=seed,kind=kind))
    if include_defaults:
        for space in ['PCA30','UMAP2']:
            for method in ['SNN','HDBSCAN_R']:add(space,method,kind='default_parity')
        return result
    for space in ['UMAP10','UMAP30','RNA_noDR']:
        for method in ['SNN','HDBSCAN_R']:add(space,method,kind='dimension')
    for space in ['PCA30','UMAP2']:
        for minpts in [25,100]:add(space,'HDBSCAN_R',minpts=minpts)
        for res in [.25,1.]:add(space,'SNN',res=res)
        for seed in [0,1,2,3]:
            for method in ['SNN','HDBSCAN_R']:add(space,method,seed=seed,kind='embedding_seed')
    assert len(result)==30
    return result
