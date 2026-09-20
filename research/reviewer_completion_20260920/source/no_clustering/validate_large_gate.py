"""Completed large-sample A2 integrity gate before an operational throttle increase."""
from pathlib import Path
from datetime import datetime, timezone
import json
import os
import subprocess
import sys

ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/no_clustering'
sys.path.insert(0,str(OUT/'source_v2'))
from a2_common import verify_source,verify_terminal,checked,sha,require_slurm,write_json


def main(sample):
    require_slurm()
    import numpy as np
    import pandas as pd
    assert sample in ['NL090','SN040']
    protocol=verify_source();results=[];seeds=[]
    tasks=json.loads((OUT/'full_tasks.json').read_text())
    index=next(i for i,t in enumerate(tasks) if t['sample']==sample)
    job='7445870_'+str(index)
    accounting=subprocess.run(['sacct','-n','-X','-P','-j',job,'--format=JobID,State,ExitCode'],capture_output=True,text=True,check=True).stdout
    records=[line.split('|') for line in accounting.splitlines() if line.strip()]
    assert any(row[0]==job and row[1]=='COMPLETED' and row[2]=='0:0' for row in records)
    for budget in protocol['budgets']:
        unit=OUT/'GBM'/sample/budget
        assert checked(unit,'fit_manifest.json','FIT_COMPLETE') and checked(unit/'evaluation')
        cfg=json.loads((unit/'config.json').read_text())
        fit=json.loads((unit/'fit_manifest.json').read_text())
        assert fit['input_signature']==cfg['input_signature']
        em=json.loads((unit/'evaluation/manifest.json').read_text())
        assert em['source_bundle_sha256']==protocol['source_bundle_sha256']
        for name,digest in em['outputs'].items():assert sha(unit/'evaluation'/name)==digest
        assert em['n_candidates']==5 and em['n_metric_rows']==15
        n_cells=protocol['inputs'][sample+'/'+budget]['n_cells']
        assert em['n_cells']==n_cells and n_cells>=10000
        seed=pd.read_csv(unit/'cellwise_seed/initial_calls.csv.gz',dtype=str,keep_default_na=False)
        assert len(seed)==n_cells and not seed.cell_id.duplicated().any()
        seeds.append(seed)
        states=[]
        for arm in protocol['arms']:
            aid=arm['id'];tm=verify_terminal(unit/'cellwise_seed',aid)
            assert fit['terminal_manifests'][aid]==sha(unit/'cellwise_seed/terminal'/aid/'terminal_manifest.json')
            with np.load(unit/'cellwise_seed/terminal'/aid/'terminal.npz') as z:
                assert np.array_equal(z['initial'],seed[aid])
                known=z['initial']!='Undecided'
                for key in ['final090','final070']:
                    assert len(z[key])==n_cells and np.array_equal(z[key][known],z['initial'][known])
            states.append({'arm':aid,'dl_status':tm['dl_status'],'training_executed':tm['training_executed']})
        metrics=pd.read_csv(unit/'evaluation/metrics.csv')
        assert len(metrics)==15 and metrics.n_cells.eq(n_cells).all()
        assert set(metrics.lambda_value)=={0,.5,1,1.5,2}
        results.append(dict(budget=budget,n_cells=n_cells,all_five_candidates_present=True,
                            terminal_and_output_hashes_valid=True,execution_states=states,
                            evaluation_manifest_sha256=sha(unit/'evaluation/manifest.json')))
    assert seeds[0].equals(seeds[1])
    report=dict(status='passed',sample=sample,job_validated=job,successful_no_OOM_exit=True,
                all_RNA_seed_identity_across_budgets=True,budgets=results,
                source_bundle_sha256=protocol['source_bundle_sha256'],validation_source_sha256=sha(Path(__file__)),
                scope='Large-sample integrity and resource gate only; not whole A2 or work-package A completion',
                job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),completed_at=datetime.now(timezone.utc).isoformat())
    write_json(OUT/(sample+'_large_validation.json'),report)
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main(sys.argv[1])
