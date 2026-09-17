#!/usr/bin/env python3
"""Capture SLURM accounting for this experiment's explicitly identified jobs only."""
import io
import json
import subprocess
from common import OUT, require_slurm, utc, write_json, sha


def collect():
    require_slurm()
    import pandas as pd
    out = OUT / 'resources'
    out.mkdir(exist_ok=True)
    groups = {
        'preparation': ['7327791','7327813','7327851'],
        'fit_pilot': ['7327993','7327995'],
        'full_fitting': ['7328008','7328106','7328647','7328175','7328934'],
        'numerical_diagnostic': ['7328984','7329127','7329541','7329589'],
        'validation_pilot': ['7327982','7327996','7328815','7329399'],
        'secondary_pilot': ['7329241','7329354'],
        'reporting_pilot': ['7328002','7328814','7328691','7329109'],
        'controller': ['7328690','7329437','7329633'],
    }
    dispatch = json.loads((OUT / 'dispatch_ledger.json').read_text())
    for stage in ['evaluation','quality']:
        groups[stage] = [v['job'] for v in dispatch[stage].values()]
    secondary = json.loads((OUT / 'singleton_dispatch_ledger.json').read_text())
    groups['secondary_singleton'] = [v['job'] for v in secondary['jobs'].values()]
    final = json.loads((OUT / 'finalization_ledger.json').read_text())
    for stage, row in final['jobs'].items():
        groups[f'delivery_{stage}'] = [row['job']] + row.get('prior_jobs', [])
    lookup = {str(job): stage for stage, jobs in groups.items() for job in jobs}
    ids = sorted(lookup)
    fields = ['JobIDRaw','JobID','JobName','State','ExitCode','ElapsedRaw','TotalCPU',
              'AllocCPUS','ReqMem','MaxRSS','MaxVMSize','Start','End','NodeList']
    records = []
    for start in range(0,len(ids),150):
        cmd = ['sacct','-n','-P','-S','2026-09-16T00:00:00','-j',','.join(ids[start:start+150]),
               '--format=' + ','.join(name + ('%100' if name in ['JobID','JobName','NodeList'] else '') for name in fields)]
        result = subprocess.run(cmd,capture_output=True,text=True,check=True)
        if result.stdout.strip():
            records.append(pd.read_csv(io.StringIO(result.stdout),sep='|',header=None,names=fields,dtype=str,keep_default_na=False))
    data = pd.concat(records,ignore_index=True).drop_duplicates('JobIDRaw')
    data['stage'] = data.JobID.str.split('_').str[0].str.split('.').str[0].map(lookup)
    assert data.stage.notna().all()
    data.to_csv(out / 'slurm_job_and_step_accounting.csv',index=False)
    allocations = data[~data.JobIDRaw.str.contains(r'\.')].copy()
    allocations['elapsed_seconds'] = pd.to_numeric(allocations.ElapsedRaw)
    allocations['allocated_cpu_hours'] = allocations.elapsed_seconds * pd.to_numeric(allocations.AllocCPUS) / 3600
    def kib(value):
        if not value:
            return float('nan')
        scale = {'K':1,'M':1024,'G':1024**2,'T':1024**3}
        return float(value[:-1])*scale[value[-1]] if value[-1] in scale else float(value)/1024
    steps = data[data.JobIDRaw.str.endswith('.batch')].copy()
    steps['allocation_raw_id'] = steps.JobIDRaw.str.removesuffix('.batch')
    steps['peak_rss_GiB'] = steps.MaxRSS.map(kib) / 1024**2
    allocations = allocations.merge(steps[['allocation_raw_id','peak_rss_GiB']],left_on='JobIDRaw',right_on='allocation_raw_id',how='left',validate='one_to_one')
    allocations.to_csv(out / 'slurm_allocation_resources.csv',index=False)
    summary = allocations.groupby('stage').agg(n_allocations=('JobIDRaw','size'),
        allocated_cpu_hours=('allocated_cpu_hours','sum'),total_wall_hours=('elapsed_seconds',lambda x:x.sum()/3600),
        largest_observed_peak_rss_GiB=('peak_rss_GiB','max'))
    summary.to_csv(out / 'slurm_resource_summary.csv')
    write_json(out / 'manifest.json',dict(timestamp=utc(),status='completed',job_groups=groups,
        note='Accounting snapshot of explicit run IDs. Reserved CPU hours are not measured CPU utilization. '
             'Active delivery/controller allocations are incomplete. FAILED fit allocations include verified singleton-scoring limitations; use scientific availability tables to interpret failures. '
             'SLURM MaxRSS is sampled and may miss transient peaks. Any unlisted exploratory work is outside this accounting scope.',
        outputs={p.name:sha(p) for p in out.glob('*.csv')}))
    return summary


if __name__ == '__main__':
    print(collect().to_string())
