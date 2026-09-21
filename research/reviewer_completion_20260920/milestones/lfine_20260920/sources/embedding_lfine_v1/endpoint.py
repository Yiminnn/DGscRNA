"""Exact frozen legacy endpoint functions, AST-extracted without module side effects."""
import ast
import common as c
def initialize():
    c.require_slurm()
    import numpy as np
    import pandas as pd
    def extract(path,names,namespace):
        tree=ast.parse(path.read_text())
        nodes=[n for n in tree.body if (isinstance(n,ast.FunctionDef) and n.name in names) or
            (isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id in names for t in n.targets))]
        assert len(nodes)==len(names)
        exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),namespace)
        return namespace
    helpers=extract(c.REFERENCE/'grid_sets.py',{'compose_lfine','lfine_targets','macro_f1_lfine','LFINE_PREFIX'},{'np':np,'pd':pd})
    exact=extract(c.REFERENCE/'v5_final_annotations.py',{'label_scope','metrics'},{'np':np,'pd':pd})
    mapping=pd.read_csv(c.REFERENCE/'panel_L1_mapping.csv',dtype=str,keep_default_na=False)
    frame=mapping[mapping.library.eq(c.FIXED)];assert len(frame) and not frame.panel.duplicated().any()
    lookup=dict(zip(frame.panel,frame.L1));helpers['panel_labels']={c.FIXED:set(lookup.values())}
    prefixes={k:list(v) for k,v in helpers['LFINE_PREFIX'].items()}
    assert prefixes==c.js(c.REFERENCE/'lfine_target_prefixes.json')
    return helpers,exact,lookup

def score(native,truth,scope,helpers,exact,lookup):
    import numpy as np
    mapped=np.asarray([lookup.get(v,'Unknown' if v in c.ABSTAIN else 'UNMAPPABLE') for v in native],dtype=object)
    values=exact['metrics'](mapped,truth,scope,helpers,c.FIXED)
    assert values['n_reference_lfine_disagreements']==0
    lf,classes,targets=scope
    # Publish counts for every observed reference class, including classes that
    # are not in the macro average; the reference does not choose predictions.
    ok=np.asarray([g in targets.get(q,()) for g,q in zip(lf,mapped)])
    perclass=[]
    for label in sorted(set(lf)):
        gold=lf==label;claimed=np.asarray([label in targets.get(q,()) for q in mapped])
        tp=int((gold&ok).sum());fn=int((gold&~ok).sum());fp=int((~ok&claimed&~gold).sum())
        perclass.append(dict(lfine_class=label,support=int(gold.sum()),included_in_macro=label in classes,
            TP=tp,FN=fn,FP=fp,F1=2*tp/(2*tp+fn+fp) if 2*tp+fn+fp else 0.))
    return {k:values[k] for k in c.FIELDS},perclass
