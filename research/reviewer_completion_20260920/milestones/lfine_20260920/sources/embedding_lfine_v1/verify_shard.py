"""Independent per-shard Lfine metric verification, saved labels only, SLURM only."""
import argparse,os
from pathlib import Path
import common as c
def run(index,space):
    c.require_slurm();contract=c.contract();c.archive_gate(contract)
    import numpy as np
    import pandas as pd
    from independent import compute
    task=c.js(c.OUT/'shards.json')[index];sample,budget=task['sample'],task['budget']
    assert space in task['spaces']
    source=c.OUT/'evaluation'/sample/budget/space;m=c.checked_outputs(source)
    out=c.OUT/'verification'/sample/budget/space
    if c.checked(out):return c.checked_outputs(out)
    out.mkdir(parents=True,exist_ok=True)
    truthpath=c.NATIVE/'evaluation_inputs'/sample/'truth.csv.gz'
    truth=pd.read_csv(truthpath,dtype=str,keep_default_na=False)
    input_info=c.js(c.NATIVE/'inputs'/sample/'input_manifest.json');assert c.sha(truthpath)==input_info['evaluation_files']['truth.csv.gz']
    mapping=pd.read_csv(c.REFERENCE/'panel_L1_mapping.csv',dtype=str,keep_default_na=False);mapping=mapping[mapping.library.eq(c.FIXED)]
    lookup=dict(zip(mapping.panel,mapping.L1));prefixes=c.js(c.REFERENCE/'lfine_target_prefixes.json')
    data=pd.read_csv(source/'metrics.csv.gz',dtype={'cutoff':str});anchors=pd.read_csv(source/'anchor_metrics.csv',dtype={'cutoff':str})
    perclass=pd.read_csv(source/'per_class.csv.gz');combined=pd.concat([data,anchors],ignore_index=True)
    assert len(data)==39 and len(anchors)==3 and not combined.duplicated(['space','route','stage']).any()
    assert data.space.eq(space).all() and set(data.stage)==set(c.STAGES)
    assert combined.terminal_valid.dtype==bool and combined.primary.dtype==bool
    assert combined.primary.eq(input_info['primary']).all() and combined.patient.eq(input_info['patient']).all()
    assert combined.n_cells.eq(len(truth)).all() and truth.cell_id.is_unique
    rows=[];total_classes=0;cache={};maximum=0.
    for row in combined.itertuples(index=False):
        if not row.terminal_valid:
            assert row.status in ['invalid','unavailable']
            for key in c.MEASURES:assert pd.isna(getattr(row,key))
            rows.append(dict(space=row.space,route=row.route,stage=row.stage,status='NA_preserved',n_cells=len(truth)));continue
        assert row.status=='completed'
        path=Path(row.predictions_path)
        if str(path) not in cache:
            assert c.sha(path)==row.predictions_sha256;cache[str(path)]=pd.read_csv(path,dtype=str,keep_default_na=False)
        pred=cache[str(path)];assert np.array_equal(pred.cell_id,truth.cell_id)
        values,counts=compute(pred[c.STAGES[row.stage]],truth,lookup,prefixes)
        for key,value in values.items():
            actual=getattr(row,key);assert np.isclose(value,actual,rtol=0,atol=1e-12,equal_nan=True),(row.space,row.route,row.stage,key,value,actual)
            if np.isfinite(value):maximum=max(maximum,abs(value-actual))
        got=perclass[(perclass.space==row.space)&(perclass.route==row.route)&(perclass.stage==row.stage)].set_index('lfine_class')
        assert not got.index.duplicated().any() and set(got.index)=={x['lfine_class'] for x in counts}
        for item in counts:
            for key in ['support','included_in_macro','TP','FN','FP','F1']:
                assert np.isclose(got.loc[item['lfine_class'],key],item[key],rtol=0,atol=1e-12)
        assert sum(x['support'] for x in counts)==len(truth)
        total_classes+=len(counts)
        rows.append(dict(space=row.space,route=row.route,stage=row.stage,status='passed',n_cells=len(truth),n_classes=len(counts),n_macro_classes=values['lfine_n_classes']))
    assert total_classes==len(perclass)
    pd.DataFrame(rows).to_csv(out/'checks.csv',index=False)
    # Parity with compact v5 native-R endpoint, at the exact same saved anchor.
    compact=c.ROOT/'results/hvg_ptc_20260916_v1/lfine_compact_20260920/sample_metrics'/f'{sample}.csv'
    old=pd.read_csv(compact,dtype={'cutoff':str});old=old[old.budget.eq(budget)&old.route.eq('UMAP2_HDBSCAN_R')&old.library.eq(c.FIXED)&old.stage.eq('terminal090')]
    assert len(old)==1
    current=anchors[anchors.stage.eq('terminal090')].iloc[0]
    assert current.predictions_sha256==old.iloc[0].predictions_sha256
    assert current.truth_sha256==old.iloc[0].truth_sha256
    for key in c.FIELDS:assert np.isclose(current[key],old.iloc[0][key],rtol=0,atol=1e-12,equal_nan=True),key
    proof=dict(status='passed',sample=sample,budget=budget,space=space,n_A1_rows=39,n_anchor_rows=3,n_per_class_rows=total_classes,
        every_observed_class_preserved=True,all_cells_retained=True,macro_support_threshold=20,
        independent_contingency_formulas=True,compact_native_anchor_parity=True,max_absolute_metric_difference=maximum,
        n_valid=int(combined.terminal_valid.sum()),n_unavailable=int((~combined.terminal_valid).sum()),
        contract_sha256=c.sha(c.OUT/'contract.json'),inputs={str(source/'manifest.json'):c.sha(source/'manifest.json'),str(compact):c.sha(compact)},
        outputs={'checks.csv':c.sha(out/'checks.csv')},job=os.environ['SLURM_JOB_ID'],completed_at=c.utc(),no_fitting=True)
    c.write(out/'manifest.json',proof);c.complete(out);return proof
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,default=int(os.environ.get('SLURM_ARRAY_TASK_ID','-1')));p.add_argument('--space');a=p.parse_args()
    assert 0<=a.index<242
    for space in ([a.space] if a.space else c.js(c.OUT/'shards.json')[a.index]['spaces']):print(run(a.index,space),flush=True)
