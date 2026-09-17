"""Export fixed curator-retained cells to sparse R inputs without exposing GT to fitting."""
import hashlib,json,os,re
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'
KEYS=['brain_GBM','breast_TNBC','colorectal','kidney_ccRCC','blood_DLBCL',
      'baron_human','muraro','segerstolpe','xin','immune_ALL_human','HCL']

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    key=KEYS[int(os.environ['SLURM_ARRAY_TASK_ID'])]
    source=json.loads((OUT/'dataset_inventory.json').read_text())[key]['path']
    a=ad.read_h5ad(source,backed='r')
    parent=OUT/'inputs'/key;parent.mkdir(parents=True,exist_ok=True)
    def is_counts(m):
        z=m[:min(500,m.shape[0]),:]
        vals=z.data if sp.issparse(z) else np.asarray(z).ravel()
        return bool(len(vals) and np.isfinite(vals).all() and vals.min()>=0 and
                    np.allclose(vals,np.rint(vals),rtol=0,atol=1e-5))
    if is_counts(a.X):
        matrix=a.X;var=a.var;layer='X';semantics='integer_counts'
    elif a.raw is not None and is_counts(a.raw.X):
        matrix=a.raw.X;var=a.raw.var;layer='raw.X';semantics='integer_counts'
    elif key=='xin':
        matrix=a.X;var=a.var;layer='X';semantics='published_RPKM_no_raw_counts'
    else:
        raise RuntimeError(f'{key}: no verified count input; inspect source before fitting')
    if key=='immune_ALL_human':
        semantics='published_counts_layer_mixed_UMI_and_full_length_quantification'
    symbol_col=next((c for c in ['feature_name','gene_symbols','gene_symbol','Symbol'] if c in var),None)
    symbols=np.array(var[symbol_col].astype(str) if symbol_col else var.index.astype(str))
    keep=np.array([bool(s.strip()) and s.lower() not in {'nan','none'} for s in symbols])
    valid=symbols[keep]
    genes=list(dict.fromkeys(valid))
    lookup={g:i for i,g in enumerate(genes)}
    merger=sp.csr_matrix((np.ones(len(valid)),(np.arange(len(valid)),[lookup[g] for g in valid])),shape=(len(valid),len(genes)))
    batch_col=next(c for c in ['donor','donor_id','Patient','sample_id','Sample.name'] if c in a.obs)
    gt_col='GT' if 'GT' in a.obs else 'cell_type'
    if key=='HCL':
        scopes=[(str(t),np.flatnonzero(a.obs.tissue_sample.astype(str).to_numpy()==t)) for t in sorted(a.obs.tissue_sample.astype(str).unique())]
    else:scopes=[('whole',np.arange(a.n_obs))]
    units=[]
    for scope,indices in scopes:
        unit=key if scope=='whole' else f'HCL__{re.sub("[^A-Za-z0-9]+","_",scope)}'
        dest=parent/unit;dest.mkdir(exist_ok=True)
        if (dest/'EXPORT_COMPLETE').exists():
            units.append(json.loads((dest/'input_manifest.json').read_text()));continue
        x=sp.csr_matrix(matrix[indices,:],dtype=np.float64)
        x=(x[:,keep]@merger).tocsr();x.eliminate_zeros();x.sort_indices()
        assert np.isfinite(x.data).all() and x.data.min()>=0
        if semantics=='integer_counts':assert np.allclose(x.data,np.rint(x.data),rtol=0,atol=1e-5)
        assert len(x.data)<2**31 and x.shape[0]<2**31
        # CSR cell x gene arrays are CSC gene x cell arrays with the same storage.
        x.data.astype('<f8').tofile(dest/'x.bin')
        x.indices.astype('<i4').tofile(dest/'i.bin')
        x.indptr.astype('<i4').tofile(dest/'p.bin')
        pd.DataFrame({'gene':genes}).to_csv(dest/'genes.csv',index=False)
        obs=a.obs.iloc[indices]
        fit=pd.DataFrame({'cell_id':obs.index.astype(str),'batch':obs[batch_col].astype(str).to_numpy()})
        assert fit.cell_id.is_unique
        fit.to_csv(dest/'cells_fit.csv',index=False)
        truth=pd.DataFrame({'cell_id':obs.index.astype(str),'truth':obs[gt_col].astype(str).to_numpy(),
                            'donor':obs[batch_col].astype(str).to_numpy()})
        for c in ['cell_type_ontology_term_id','tissue','tissue_sample','stage','author_cell_type','Sample ID']:
            if c in obs:truth[c]=obs[c].astype(str).to_numpy()
        truth.to_csv(dest/'evaluation_only.csv.gz',index=False)
        record=dict(unit=unit,dataset=key,scope=scope,path=str(dest),source=source,
            expression_layer=layer,input_semantics=semantics,n_cells=x.shape[0],n_genes=x.shape[1],nnz=x.nnz,
            duplicate_symbols_aggregated=int(len(valid)-len(genes)),batch_column=batch_col,
            batch_sizes=fit.batch.value_counts().to_dict(),curator_cell_set_preserved=True,
            evaluation_labels_excluded_from_fit=True,job=os.environ['SLURM_JOB_ID'],
            script_sha256=sha(__file__),files={n:sha(dest/n) for n in ['x.bin','i.bin','p.bin','genes.csv','cells_fit.csv']})
        (dest/'input_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
        (dest/'EXPORT_COMPLETE').write_text(sha(dest/'input_manifest.json')+'\n')
        units.append(record)
        print(unit,x.shape,'batches',record['batch_sizes'],semantics,flush=True)
    (parent/'units.json').write_text(json.dumps(units,indent=2)+'\n')
    a.file.close()

if __name__=='__main__':run()
