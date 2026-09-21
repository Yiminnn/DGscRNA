"""Copied C ordered reductions; no original imports or execution."""
SPEC=['method','budget','route','library','cutoff']

MEASURES=['lfine_macroF1','coverage','legacy_called_coverage','abstain_rate',
    'offvocab_rate','acc_on_called','lfine_scored_class_cell_fraction','marker_vocab_oracle_upper_bound']

def ordered_patient_means(data,folds):
    """One canonical IEEE-754 reduction order, used for exact lexical tie ranking."""
    ordered=data.sort_values(['sample']+SPEC,kind='stable')
    patient=ordered.groupby(['patient']+SPEC,sort=True)[MEASURES].mean().reset_index()
    return patient.merge(folds,on='patient',validate='many_to_one').sort_values(['patient']+SPEC,kind='stable').reset_index(drop=True)

def ordered_training_rank(patient,fold):
    training=patient[patient.fold.ne(fold)].sort_values(['patient']+SPEC,kind='stable')
    scores=training.groupby(SPEC,sort=True)[MEASURES].mean().reset_index()
    return scores.sort_values(['lfine_macroF1']+SPEC,ascending=[False]+[True]*len(SPEC),kind='stable').reset_index(drop=True)
