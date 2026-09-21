"""Synthetic endpoint boundary checks; only execute inside an authorized SLURM job."""
import common as c
def run():
    c.require_slurm();c.contract()
    import numpy as np
    import pandas as pd
    from endpoint import initialize,score
    from independent import compute
    helpers,exact,_=initialize()
    labels=['Malignant_OPC']*20+['OPC']*20+['TAM_MD']*20+['TAM_MG']*19+['Other']*2
    truth=pd.DataFrame(dict(L1=['Malignant' if x.startswith('Malignant') else 'TAM' if x.startswith('TAM') else x for x in labels],
        L3=['' if x.startswith('Malignant') else x for x in labels],MalState=['OPC' if x.startswith('Malignant') else '' for x in labels],
        lfine_original=labels,Lfine=labels))
    lookup={'panel_OPC':'OPC','panel_TAM':'TAM','panel_Other':'Other','panel_wrong':'UNMAPPABLE'}
    helpers['panel_labels']={c.FIXED:set(lookup.values())}
    prefixes=c.js(c.REFERENCE/'lfine_target_prefixes.json')
    scope=exact['label_scope'](truth,helpers)
    assert set(scope[1])=={'Malignant_OPC','OPC','TAM_MD'} and 'TAM_MG' in scope[2]['TAM']
    cases={
        'compatible_all_correct':['panel_OPC']*40+['panel_TAM']*39+['panel_Other']*2,
        'all_Unknown':['Unknown']*81,'all_Undecided':['Undecided']*81,
        'off_vocabulary':['panel_wrong']*81,'wrong_broad_claim':['panel_TAM']*81,
        'mixed_abstention':['Noise','None','nan','']*20+['unknown_native_panel']}
    records=[]
    for name,native in cases.items():
        values,pc=score(np.asarray(native),truth,scope,helpers,exact,lookup);got,counts=compute(native,truth,lookup,prefixes)
        for key in c.FIELDS:assert np.isclose(values[key],got[key],rtol=0,atol=1e-12,equal_nan=True),(name,key)
        assert pc==counts
        if name=='compatible_all_correct':assert values['lfine_macroF1']==1
        if name in ['all_Unknown','all_Undecided','off_vocabulary']:assert values['lfine_macroF1']==0
        records.append(dict(case=name,passed=True))
    small=truth.iloc[:10].copy();small_scope=exact['label_scope'](small,helpers);assert not small_scope[1]
    values,_=score(np.asarray(['panel_OPC']*10),small,small_scope,helpers,exact,lookup);got,_=compute(['panel_OPC']*10,small,lookup,prefixes)
    assert np.isnan(values['lfine_macroF1']) and np.isnan(got['lfine_macroF1']) and values['coverage']==1
    records.append(dict(case='no_support20_class_is_NA_not_zero',passed=True))
    return dict(status='passed',n_cases=len(records),cases=records)
if __name__=='__main__':print(run())
