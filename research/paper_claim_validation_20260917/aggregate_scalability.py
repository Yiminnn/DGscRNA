"""Actual elapsed time and scheduler peak RSS, including failed attempts."""
import json
import os
import re
from common import OUT, require_slurm, checked, complete, sha, write_json, utc

def run():
    require_slurm()
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scalability_dispatch import SIZES,PILOTS,dest,accounting
    state=json.loads((OUT/'scalability_dispatch_state.json').read_text())
    assert state['completed']==15 and not state['review_required']
    target=OUT/'scalability_summary';target.mkdir(exist_ok=True)
    rows=[];accounting_rows=[];sources={}
    for method in PILOTS:
        for n in SIZES:
            source=dest(method,n);assert checked(source)
            m=json.loads((source/'manifest.json').read_text());sources[f'{method}/{n}']=sha(source/'manifest.json')
            acc=accounting(m['job']);peak=[]
            for r in acc:
                accounting_rows.append(dict(method=method,n_cells=n,**r))
                if r.get('MaxRSS'):
                    value=re.fullmatch(r'([0-9.]+)([KMGTP]?)',r['MaxRSS']);assert value
                    scale={'':1,'K':1024,'M':1024**2,'G':1024**3,'T':1024**4,'P':1024**5}[value[2]]
                    peak.append(float(value[1])*scale/1024**3)
            assert peak,(method,n,'No actual MaxRSS accounting')
            assert m.get('cold_preparation',m.get('cold_input',False)),(method,n,'Not a cold run')
            rows.append(dict(method=method,n_cells=n,pipeline_seconds=m['pipeline_elapsed_seconds'],
                job_step_peak_RSS_GiB=max(peak),hardware=m['hardware'],CPUs=m['SLURM_CPUS_PER_TASK'],job=m['job'],
                scope='Fixed marker default' if method!='scDeepSort' else 'Published pretrained Brain GNN'))
    frame=pd.DataFrame(rows);frame.to_csv(target/'cold_pipeline_resource_curve.csv',index=False)
    pd.DataFrame(accounting_rows).to_csv(target/'slurm_job_and_step_accounting.csv',index=False)
    fig,axs=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for method,color in zip(PILOTS,['#0072B2','#D55E00','#009E73']):
        data=frame[frame.method==method].sort_values('n_cells')
        axs[0].plot(data.n_cells/1000,data.pipeline_seconds/60,'o-',color=color,label=method)
        axs[1].plot(data.n_cells/1000,data.job_step_peak_RSS_GiB,'o-',color=color,label=method)
    for ax in axs:ax.set_xlabel('Cells in one workflow run (thousands)');ax.legend(frameon=False)
    axs[0].set_ylabel('Cold pipeline wall time (minutes)');axs[1].set_ylabel('Recorded job/step peak RSS (GiB)')
    fig.suptitle('Nested pooled-count resource test: 10k to 120k cells\nCPU allocation recorded; distinct marker/reference information; no biological batch inference')
    for ext in ['png','pdf']:fig.savefig(target/f'resource_curves.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)
    note='Each point starts from counts without a completed-result cache. Operating-system file-cache state is uncontrolled; cold input does not mean a cold machine. Each point is one CPU-only run, not a timing replicate or a sum of independent samples; GPU memory is not applicable. The same nested pooled-cell inputs are used. DG-scRNA runs the original HVG2000/PCA30/UMAP2/HDBSCAN50 workflow and fixed glioma-marker terminal DL; its shared scoring implementation also calculates the other marker-density tables. SCINA receives the same fixed glioma marker library and needs no clustering. scDeepSort uses its published human Brain pretrained model. Each submission requested four task CPUs. Cluster memory policy can increase the actual CPU allocation; the table reports those actual allocations and hardware, while computational thread caps remain method-specific. SCINA includes normalization and scDeepSort includes input conversion. DG plotting is excluded from pipeline wall time. MaxRSS is the scheduler-reported job/step measurement, not a claim about summed simultaneous RSS of all forked children. The run does not establish that pooling patients removes batch effects, nor does a completed point alone establish resource superiority.\n'
    (target/'RESOURCE_INTERPRETATION.md').write_text(note)
    write_json(target/'manifest.json',dict(status='completed',n_cold_runs=15,largest_single_run=120000,
        timing_replicates=1,source_manifests=sources,files={p.name:sha(p) for p in target.iterdir() if p.suffix in ['.csv','.png','.pdf','.md']},
        job=os.environ['SLURM_JOB_ID'],source_sha256=sha(__file__),completed_at=utc()))
    complete(target)

if __name__=='__main__':run()
