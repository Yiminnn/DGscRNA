#!/usr/bin/env python3
"""Independent saved-call arithmetic and full roster verification, no v5 import."""
from collections import Counter
from pathlib import Path
import json
import os
from evaluate import ROOT,HERE,OLD,NATIVE,COMPACT,OUT,LIBRARIES,PILOTS,sha,write_json,utc,verify_freeze


def main():
    assert os.environ.get('SLURM_JOB_ID')
    verify_freeze()
    import numpy as np
    import pandas as pd
    dest=OUT/'independent_verification'
    dest.mkdir(exist_ok=False)
    full=OUT/'full'
    manifest=json.loads((full/'validation.json').read_text())
    assert manifest['status']=='passed' and manifest['n_rows']==408
    for name,digest in manifest['files'].items(): assert sha(full/name)==digest,name
    artifacts=json.loads((full/'artifact_hashes.json').read_text())
    for name,digest in artifacts.items(): assert sha(ROOT/name)==digest,name
    result=pd.read_csv(full/'terminal_metrics.csv.gz')
    perclass=pd.read_csv(full/'per_class_counts.csv.gz')
    keys=['task','sample','budget','control','library','model_seed','epochs','stage']
    expected=[]
    for cfg in json.loads((OLD/'tasks.json').read_text()):
        if cfg['task']=='neighbors':
            for library in LIBRARIES:
                for stage in ['terminal090','terminal070']:
                    expected.append(('neighbors',cfg['sample'],cfg['budget'],cfg['name'],library,42,10,stage))
        else:
            for epoch in [5,10,20,30]:
                for stage in ['terminal090','terminal070']:
                    expected.append(('learning',cfg['sample'],cfg['budget'],'learning',LIBRARIES[0],cfg['model_seed'],epoch,stage))
    got=list(result[keys].itertuples(index=False,name=None))
    assert len(expected)==len(set(expected))==len(got)==len(set(got))==408
    assert set(expected)==set(got)
    assert Counter(zip(result.task,result.stage))==Counter({('neighbors','terminal090'):132,('neighbors','terminal070'):132,('learning','terminal090'):72,('learning','terminal070'):72})
    assert not any(c in result for c in ['macroF1_present','accuracy','malignant_F1'])
    assert result.terminal_valid.all() and not result.backup.any()
    prefix=json.loads((COMPACT/'lfine_target_prefixes.json').read_text())
    assert prefix==json.loads((full/'lfine_target_prefixes.json').read_text())
    mapping=pd.read_csv(NATIVE/'markers/panel_L1_mapping.csv',dtype=str,keep_default_na=False)
    maps={library:dict(zip(group.panel,group.L1)) for library,group in mapping.groupby('library')}
    truth={s:pd.read_csv(NATIVE/'evaluation_inputs'/s/'truth.csv.gz',dtype=str,keep_default_na=False) for s in PILOTS}
    raw_abstain={'Unknown','Undecided','Noise','nan','None',''}
    proof=[]
    for row in result.itertuples(index=False):
        t=truth[row.sample]
        lf=t.lfine_original.to_list()
        count=Counter(lf)
        classes=[c for c,n in count.items() if n>=20 and c not in {'Other','nan'}]
        targets={parent:{c for c in count if any(c.startswith(p) for p in starts)} for parent,starts in prefix.items()}
        pred=pd.read_csv(ROOT/row.predictions_path,dtype=str,keep_default_na=False)
        assert list(pred.cell_id)==list(t.cell_id)
        raw=pred['final'+row.stage[-3:]].to_list()
        mapped=[maps[row.library].get(call,'Unknown' if call in raw_abstain else 'UNMAPPABLE') for call in raw]
        # Count confusion contributions in a single cell loop, independently of
        # the vectorized v5/evaluator implementation and its macro helper.
        tp=Counter();fn=Counter();fp=Counter();hits=set();correct_called=0
        abstain_count=offvocab_count=called_count=0;distinct=set()
        for label,call in zip(lf,mapped):
            compatible=targets.get(call,set())
            correct=label in compatible
            abstain=call in raw_abstain
            offvocab=not abstain and call not in prefix
            called=not abstain and not offvocab
            abstain_count+=int(abstain);offvocab_count+=int(offvocab);called_count+=int(called)
            if called:
                distinct.add(call);correct_called+=int(correct)
            if correct:
                tp[label]+=1;hits.add(label)
            else:
                fn[label]+=1
                for candidate in compatible:
                    if candidate!=label: fp[candidate]+=1
        fs=[]
        pc=perclass[perclass.row_key.eq(row.row_key)].set_index('lfine_class')
        assert set(pc.index)==set(classes) and pc.index.is_unique
        for label in classes:
            denominator=2*tp[label]+fp[label]+fn[label]
            f=2*tp[label]/denominator if denominator else 0.0
            fs.append(f)
            r=pc.loc[label]
            assert (int(r.TP),int(r.FN),int(r.FP),int(r.support))==(tp[label],fn[label],fp[label],count[label])
            assert abs(r.F1-f)<1e-12
        reachable=set().union(*(targets.get(parent,set()) for parent in set(maps[row.library].values())))
        bad={'nan','','None','NA'}
        composed=['Malignant_'+m if m not in bad else l3 if l3 not in bad else l1 for l1,l3,m in zip(t.L1,t.L3,t.MalState)]
        n=len(lf)
        values=dict(lfine_macroF1=sum(fs)/len(fs) if fs else np.nan,
            coverage=(n-abstain_count)/n,legacy_called_coverage=called_count/n,
            abstain_rate=abstain_count/n,offvocab_rate=offvocab_count/n,
            acc_on_called=correct_called/called_count if called_count else np.nan,
            n_distinct_calls=len(distinct),n_classes_hit=len(hits&set(classes)),
            lfine_n_classes=len(classes),lfine_scored_class_cell_fraction=sum(count[c] for c in classes)/n,
            marker_vocab_oracle_upper_bound=len(set(classes)&reachable)/len(classes) if classes else np.nan,
            n_reference_lfine_disagreements=sum(a!=b for a,b in zip(lf,composed)))
        for field,value in values.items():
            np.testing.assert_allclose(value,getattr(row,field),rtol=0,atol=1e-12,equal_nan=True,err_msg=row.row_key+'/'+field)
        proof.append(dict(row_key=row.row_key,n_cells=n,n_fields=len(values),**values))
    pd.DataFrame(proof).to_csv(dest/'independently_recomputed_metrics.csv.gz',index=False)
    verify_freeze()
    validation=dict(status='passed',at=utc(),job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),
        n_rows=408,n_neighbor_rows=264,n_learning_rows=144,n_fields_per_row=12,
        independent_cell_loop_arithmetic=True,independent_full_roster=True,thresholds=['0.90','0.70'],
        epochs=[5,10,20,30],seeds=[0,1,42],n_samples=3,n_budgets=2,
        no_provider_metric_function_import=True,no_fitting=True,no_best_selection=True,
        no_L1_performance_values=True,source_sha256=sha(__file__),full_validation_sha256=sha(full/'validation.json'),
        protocol_sha256=sha(HERE/'PROTOCOL.md'),files={p.name:sha(p) for p in dest.iterdir() if p.is_file()})
    write_json(dest/'validation.json',validation)
    state=json.loads((OUT/'status.json').read_text())
    state.update(status='validated',updated_at=utc(),completed=408,remaining=0,
        summary='All408 saved-terminal Lfine compatible-target-set rows independently verified; three GBM pilots only.',
        validation_file=str(dest/'validation.json'))
    state['jobs'].append(dict(job=validation['job'],step=validation['step']))
    state['evidence'].append(str(dest/'validation.json'))
    write_json(OUT/'status.json',state)
    print(json.dumps(validation,indent=2),flush=True)


if __name__=='__main__':main()
