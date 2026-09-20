"""A2 fitting entrypoint: no truth labels enter cellwise seeding or refinement."""
from pathlib import Path
import argparse
import json
import os
import shutil
import subprocess
import sys

from a2_common import CODE, OUT, REFERENCE, RSCRIPT, configure, verify_source, verify_terminal, require_slurm, sha, utc, write_json, complete, checked


def finish(route, only_arm=None):
    import terminal
    terminal.OUT = OUT
    os.environ['DGSCRNA_DL_CACHE_ROOT'] = str(OUT/'DL_cache')
    terminal.finish_route(route, only_arm)


def fit(sample, budget):
    cfg = configure(sample, budget);dest = Path(cfg['dest']);route = Path(cfg['route_dir'])
    if checked(dest, 'fit_manifest.json', 'FIT_COMPLETE'):
        prior = json.loads((dest/'fit_manifest.json').read_text())
        assert prior['input_signature'] == cfg['input_signature']
        for arm in cfg['arms']:
            verify_terminal(route, arm['id'])
            assert prior['terminal_manifests'][arm['id']] == sha(route/'terminal'/arm['id']/'terminal_manifest.json')
        return cfg
    subprocess.run([RSCRIPT, str(CODE/'score_cellwise.R'), str(dest/'config.json')], check=True)
    sm = json.loads((route/'score_manifest.json').read_text())
    assert checked(route, 'score_manifest.json', 'SCORE_COMPLETE')
    assert sm['input_signature'] == cfg['input_signature']
    assert set(sm['arms']) == {arm['id'] for arm in cfg['arms']}
    finish(route)
    manifests = {}
    for arm in cfg['arms']:
        verify_terminal(route, arm['id'])
        manifests[arm['id']] = sha(route/'terminal'/arm['id']/'terminal_manifest.json')
    write_json(dest/'fit_manifest.json', dict(status='fit_complete_evaluation_pending', input_signature=cfg['input_signature'],
        sample=sample, budget=budget, seed_mechanism_replacement=True, reference_labels_used_for_fit=False,
        terminal_manifests=manifests, source_bundle_sha256=cfg['source_bundle_sha256'],
        job=os.environ['SLURM_JOB_ID'], completed_at=utc()))
    complete(dest, 'fit_manifest.json', 'FIT_COMPLETE')
    return cfg


def regression():
    import numpy as np
    import pandas as pd
    protocol = verify_source();sample='TKU4163';budget='hvg2000'
    prep=REFERENCE/'GBM'/sample/budget;source=prep/'UMAP2_HDBSCAN_R'
    dest=OUT/'regression';route=dest/'UMAP2_HDBSCAN_R';route.mkdir(parents=True,exist_ok=True)
    if checked(dest):
        assert json.loads((dest/'manifest.json').read_text())['source_bundle_sha256']==sha(CODE/'SOURCE_MANIFEST.json')
        return
    shutil.copy2(source/'clusters.csv',route/'clusters.csv')
    cfg=dict(sample=sample,budget=budget,prep=str(prep),dest=str(dest),
             frozen_protocol_sha256=sha(OUT/'protocol.json'),
             conditions=[dict(route='UMAP2_HDBSCAN_R',dest=str(route),method='HDBSCAN_R',k=None)])
    write_json(dest/'config.json',cfg)
    subprocess.run([RSCRIPT,str(CODE/'legacy/score_original_anchor.R'),str(dest/'config.json')],check=True)
    sm=json.loads((source/'score_manifest.json').read_text())
    aid=next(a for a,m in sm['arms'].items() if m['library']=='CM2_glioma_other' and m['cutoff']=='mean')
    new=pd.read_csv(route/'initial_calls.csv.gz',dtype=str,keep_default_na=False)
    old=pd.read_csv(source/'initial_calls.csv.gz',dtype=str,keep_default_na=False)
    assert np.array_equal(new.cell_id,old.cell_id)
    assert np.array_equal(new.L00_mean,old[aid]),'Original R scoring labels differ'
    finish(route,'L00_mean');tm=verify_terminal(route,'L00_mean')
    with np.load(route/'terminal/L00_mean/terminal.npz') as a,np.load(source/f'terminal/{aid}/terminal.npz') as b:
        assert set(a.files)==set(b.files)
        for key in a.files:
            if key=='probabilities':np.testing.assert_allclose(a[key],b[key],rtol=1e-6,atol=1e-7)
            elif key=='confidence_rounded':np.testing.assert_allclose(a[key],b[key],rtol=0,atol=0,equal_nan=True)
            else:assert np.array_equal(a[key],b[key]),key
    write_json(dest/'manifest.json',dict(status='passed',initial_labels_exact=True,terminal_arrays_labels_splits_exact=True,
        probability_rtol=1e-6,probability_atol=1e-7,training_executed=tm['training_executed'],dl_status=tm['dl_status'],
        original_score_manifest_sha256=sha(source/'score_manifest.json'),
        original_terminal_manifest_sha256=sha(source/f'terminal/{aid}/terminal_manifest.json'),
        source_bundle_sha256=sha(CODE/'SOURCE_MANIFEST.json'),job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest)
    print('A2_ORIGINAL_REGRESSION_PASSED',tm['dl_status'],flush=True)


def main():
    require_slurm()
    import terminal
    terminal.threads()
    parser=argparse.ArgumentParser();parser.add_argument('tasks');parser.add_argument('mode',nargs='?',default='fit')
    args=parser.parse_args();tasks=json.loads(Path(args.tasks).read_text())
    task=tasks[int(os.environ.get('SLURM_ARRAY_TASK_ID','0'))]
    try:
        if args.mode=='regression':
            regression();return
        assert checked(OUT/'regression'),'Original score and DL parity must pass before A2 fitting'
        for budget in task.get('budgets',['hvg2000','hvg5000']):
            cfg=fit(task['sample'],budget)
            subprocess.run([sys.executable,str(CODE/'evaluate.py'),cfg['dest']],check=True)
        print('A2_SAMPLE_COMPLETE',task['sample'],flush=True)
    except Exception as exc:
        write_json(OUT/'failures'/(os.environ['SLURM_JOB_ID']+'_'+os.environ.get('SLURM_ARRAY_TASK_ID','0')+'.json'),
                   dict(task=task,mode=args.mode,error=repr(exc),job=os.environ['SLURM_JOB_ID'],time=utc()))
        raise


if __name__=='__main__':main()
