"""Independent compatible-set count formulas; no legacy endpoint function calls."""
from collections import Counter
import common as c
def compute(native,truth,lookup,prefixes):
    import numpy as np
    fine=truth.lfine_original.astype(str).to_numpy()
    missing={'nan','','None','NA'}
    composed=[('Malignant_'+m) if m not in missing else (l3 if l3 not in missing else l1)
        for l1,l3,m in zip(truth.L1,truth.L3,truth.MalState)]
    assert list(fine)==composed
    counts=Counter(fine);classes=[x for x,n in counts.items() if n>=20 and x not in ['Other','nan']]
    targets={parent:{label for label in counts if any(label.startswith(prefix) for prefix in pref)} for parent,pref in prefixes.items()}
    parents=[lookup.get(label,'Unknown' if label in c.ABSTAIN else 'UNMAPPABLE') for label in native]
    contingency=Counter(zip(fine,parents));pc=[]
    for label,support in sorted(counts.items()):
        tp=sum(n for (gold,parent),n in contingency.items() if gold==label and label in targets.get(parent,set()))
        fp=sum(n for (gold,parent),n in contingency.items() if gold!=label and gold not in targets.get(parent,set()) and label in targets.get(parent,set()))
        fn=support-tp
        pc.append(dict(lfine_class=label,support=support,included_in_macro=label in classes,TP=tp,FN=fn,FP=fp,
            F1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0.))
    abstain=sum(parent in c.ABSTAIN for parent in parents)
    called=[(gold,parent) for gold,parent in zip(fine,parents) if parent not in c.ABSTAIN and parent in prefixes]
    offvocab=len(parents)-abstain-len(called)
    correct=sum(gold in targets[parent] for gold,parent in called)
    reachable=set().union(*(targets.get(parent,set()) for parent in set(lookup.values())))
    values=dict(lfine_macroF1=sum(row['F1'] for row in pc if row['included_in_macro'])/len(classes) if classes else np.nan,
        coverage=1-abstain/len(fine),legacy_called_coverage=len(called)/len(fine),abstain_rate=abstain/len(fine),offvocab_rate=offvocab/len(fine),
        acc_on_called=correct/len(called) if called else np.nan,n_distinct_calls=len({parent for _,parent in called}),
        n_classes_hit=sum(row['TP']>0 for row in pc if row['included_in_macro']),lfine_n_classes=len(classes),
        lfine_scored_class_cell_fraction=sum(counts[k] for k in classes)/len(fine),
        marker_vocab_oracle_upper_bound=sum(k in reachable for k in classes)/len(classes) if classes else np.nan,
        n_reference_lfine_disagreements=0)
    return values,pc
