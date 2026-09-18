"""One native-R workflow on one large pooled input; explicit resource-only scope."""
import json
import os
from pathlib import Path
import platform
import sys
import time
from common import OUT, require_slurm, write_json, sha, complete, checked, utc
def run(n):
    require_slurm()
    import pipeline
    n=int(n);sample=f'SCALE_{n}';dest=OUT/'scalability'/str(n);dest.mkdir(parents=True,exist_ok=True)
    if checked(dest):return
    src=OUT/'GBM'/sample/'hvg2000'
    existing=(src/'PREPARED').exists();start=time.monotonic()
    pipeline.unit(sample,'hvg2000',['UMAP2_HDBSCAN_R'],pilot=True)
    elapsed=time.monotonic()-start
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    coords=pd.read_csv(src/'UMAP2.csv',index_col=0);cl=pd.read_csv(src/'UMAP2_HDBSCAN_R/clusters.csv',dtype=str)
    assert list(coords.index)==list(cl.cell_id) and len(cl)==n
    fig,ax=plt.subplots(figsize=(7,6),layout='constrained')
    for i,label in enumerate(sorted(cl.cluster.unique())):
        take=cl.cluster.eq(label).to_numpy()
        ax.scatter(coords.iloc[take,0],coords.iloc[take,1],s=.3,c=['#bdbdbd' if label=='0' else plt.get_cmap('tab20')(i%20)],linewidths=0,rasterized=True)
    ax.set(title=f'{n:,} pooled cells: single-run resource measurement\nOriginal HVG2000 / UMAP2 / HDBSCAN50; noise in gray',xlabel='UMAP1',ylabel='UMAP2')
    for ext in ['png','pdf']:fig.savefig(dest/f'clustering.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)
    td=src/'UMAP2_HDBSCAN_R/terminal/L00_mean'
    tm=json.loads((td/'terminal_manifest.json').read_text())
    write_json(dest/'manifest.json',dict(status='completed',n_cells=n,pipeline_elapsed_seconds=elapsed,plotting_excluded_from_pipeline_timer=True,
        cold_preparation=not existing,DL_status=tm['dl_status'],cached_DL=tm.get('cache_reused'),
        annotation='Fixed CM2_glioma_other / mean -> terminal0.90; all16library density tables are produced by the shared score runner, but only the fixed-library terminal is fitted.',
        source_scope='Uniform nested pooled GBM counts, uncorrected. Runtime only; no accuracy or batch-mixing inference.',
        hardware=platform.node(),SLURM_CPUS_PER_TASK=os.environ.get('SLURM_CPUS_PER_TASK'),
        scheduler_peak_memory='Read job/step MaxRSS from completed sacct accounting; not estimated from input file size.',
        source_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest)

if __name__=='__main__':run(sys.argv[1])
