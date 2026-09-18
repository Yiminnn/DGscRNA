"""Independent checks of the R pilot before releasing the full cohort."""
import json
import os
from pathlib import Path
from common import ROOT, OUT, INPUTS, require_slurm, checked, sha, write_json, complete, utc

def run():
    require_slurm()
    import anndata as ad
    import numpy as np
    import pandas as pd
    rows=[]
    for sample in (OUT/'protocol/pilot_samples.txt').read_text().split():
        density=json.loads((OUT/'verification/pilot_density'/f'{sample}.json').read_text())
        assert density['status']=='passed' and all(v['exact'] for v in density['checks'].values())
        prep=OUT/'GBM'/sample/'hvg2000'
        assert checked(prep,'prepare_manifest.json','PREPARED')
        pm=json.loads((prep/'prepare_manifest.json').read_text())
        assert pm['features']['DL']==2000 and pm['features']['geometry']<=2000
        assert pm['features']['geometry']>=31 and not pm['reference_labels_used_for_fitting']
        assert pm['correction']=='no_correction_single_batch'
        im=json.loads((OUT/'inputs'/sample/'input_manifest.json').read_text())
        a=ad.read_h5ad(INPUTS/sample/'counts_gene_filtered.h5ad')
        cells=pd.read_csv(prep/'cells.csv').cell_id
        assert list(cells)==list(a.obs_names)
        gm=pd.read_csv(prep/'Seurat_gene_names.csv')
        assert list(gm.source)==list(a.var_names)
        assert list(gm.Seurat)==(prep/'scoring_features.txt').read_text().splitlines()
        assert pm['features']['scoring']==a.n_vars
        genes=(prep/'DL_features.txt').read_text().splitlines()
        positions=pd.Index(gm.Seurat).get_indexer(genes)
        assert (positions>=0).all()
        dense=a.X[:,positions].toarray().astype('float64')
        dense*=10000/np.asarray(a.X.sum(axis=1)).ravel()[:,None]
        np.log1p(dense,out=dense)
        dl=np.memmap(pm['DL_binary'],mode='r',dtype='<f4',shape=(a.n_obs,2000))
        np.testing.assert_allclose(dl,dense,rtol=1e-6,atol=1e-6)
        for name,width in [('PCA30.csv',30),('UMAP2.csv',2)]:
            z=pd.read_csv(prep/name,index_col=0)
            assert list(z.index)==list(cells) and z.shape==(len(cells),width)
            assert np.isfinite(z.to_numpy()).all()
        for route in ['PCA30_SNN','UMAP2_HDBSCAN_R']:
            src=prep/route;dest=src/'terminal/L00_mean'
            assert checked(src,'score_manifest.json','SCORE_COMPLETE')
            assert checked(dest,'terminal_manifest.json','TERMINAL_COMPLETE')
            tm=json.loads((dest/'terminal_manifest.json').read_text())
            assert tm['terminal_valid'] and not tm['reference_labels_used_for_fit']
            assert tm['known_labels_unchanged'] and tm['threshold_reconstructed']
            cl=pd.read_csv(src/'clusters.csv')
            assert list(cl.cell_id)==list(cells)
            rows.append(dict(sample=sample,route=route,n_cells=len(cells),DL_features=2000,
                normalized_DL_reconstructed=True,cell_order_exact=True,
                archived_density_all_libraries_exact=True,
                no_reference_labels_in_fit=True,dl_status=tm['dl_status'],terminal_valid=True))
    dest=OUT/'verification';dest.mkdir(exist_ok=True)
    pd.DataFrame(rows).to_csv(dest/'pilot_checks.csv',index=False)
    write_json(dest/'pilot_manifest.json',dict(status='passed',checks=rows,
        acceptance='Correct inputs, native stage genes, normalized DL reconstruction, complete terminal states, cell identities. No performance cutoff.',
        source_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest,'pilot_manifest.json','PILOT_PASSED')
    print('PILOT_PASSED',len(rows),flush=True)

if __name__=='__main__':run()
