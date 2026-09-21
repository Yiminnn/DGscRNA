"""Synthetic full-roster selection test; never writes the real result notebook."""
from pathlib import Path
from types import SimpleNamespace
import hashlib, importlib.util, json, os
CODE=Path(__file__).resolve().parent
ROOT=CODE.parents[2]
OUT=ROOT/'results/hvg_ptc_20260916_v1/package_reference_rerun_20260920/embedding_summary_validation/synthetic_v1'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result

def main():
 assert os.environ.get('SLURM_JOB_ID')
 assert not OUT.exists(),'Preserve synthetic validation attempts'
 import numpy as np
 import pandas as pd
 from scipy.stats import wilcoxon
 numerical=module('A1_numerical_synthetic',CODE/'numerical.py')
 wrapper=module('A1_wrapper_synthetic',CODE/'run.py')
 wrapper.source_contract()
 folds=pd.read_csv(CODE/'protocol/patient_folds.csv')
 spec=json.loads((CODE/'protocol/embedding.json').read_text());protocol=json.loads((CODE/'protocol/lfine.json').read_text())
 ineligible_patient=folds[folds.primary].patient.iloc[0]
 rows=[];anchors=[];eligibility={}
 for i,row in folds.iterrows():
  eligible=row.patient!=ineligible_patient
  eligibility[row['sample']]=dict(sample=row['sample'],patient=row.patient,primary=bool(row.primary),lfine_metric_eligible=eligible)
  for budget in spec['budgets']:
   for stage in ['marker_only','terminal070','terminal090']:
    anchor=dict(sample=row['sample'],patient=row.patient,primary=bool(row.primary),budget=budget,stage=stage,space='PCA30_UMAP2',method='HDBSCAN_R',k=0,lfine_metric_eligible=eligible)
    anchor.update({name:(.41+.001*i if name=='lfine_macroF1' and eligible else np.nan if name=='lfine_macroF1' else .8) for name in numerical.MEASURES})
    anchors.append(anchor)
    for s,space in enumerate(spec['spaces']):
     for m,method in enumerate(['KMeans','GMM','HDBSCAN_R']):
      for k in ([0] if method=='HDBSCAN_R' else spec['K']):
       value=.42+.001*i+.003*s+.004*m+.01*(budget=='hvg5000')+.008*(stage=='terminal090')
       if not(space=='noDR' and method=='KMeans'):
        value+=.025*(k==spec['K'][int(row.fold)%len(spec['K'])])+.0001*k
       candidate=dict(anchor,space=space,method=method,k=k,route=f'{space}_{method}_{k}')
       candidate.update({name:(value if eligible else np.nan) if name=='lfine_macroF1' else .8 for name in numerical.MEASURES})
       rows.append(candidate)
 data3=pd.DataFrame(rows);anchor3=pd.DataFrame(anchors)
 data=data3[data3.stage.ne('marker_only')].copy();anchor=anchor3[anchor3.stage.ne('marker_only')].copy()
 assert len(data)==44044 and len(anchor)==484
 choices,selected,patients,stats=numerical.select(data,anchor,eligibility,folds,spec,protocol)
 independent=numerical.independently_select(data,anchor,eligibility,folds,spec,protocol)
 keys=[['cohort','budget','space','method','fold'],['cohort','sample','budget','space','method','stage'],['cohort','patient','budget','space','method','stage'],['cohort','budget','space','method']]
 for left,right,key in zip([choices,selected,patients,stats],independent,keys):numerical.equal(left,right,key)
 # Execute the original scientific body unchanged with its original three stages;
 # compare its terminal endpoints to the two-endpoint adapter.
 original=(CODE/'original/select_lfine.py').read_text()
 body=original[original.index('    selected_frames=[];choices=[];patients_all=[];paired=[];rng='):original.index('    out.mkdir(parents=True,exist_ok=True)')]
 text='def original(allmetrics,anchors,eligibility,folds,spec,protocol):\n    foldmap=folds.drop_duplicates("patient").set_index("patient").fold.to_dict()\n'+body+'    return pd.DataFrame(choices),pd.concat(selected_frames),pd.concat(patients_all),stats\n'
 namespace=dict(np=np,pd=pd,wilcoxon=wilcoxon,c=SimpleNamespace(MEASURES=numerical.MEASURES))
 exec(compile(text,'<unchanged-original-A1-selection>','exec'),namespace)
 expected=namespace['original'](data3,anchor3,eligibility,folds,spec,protocol)
 for index,(left,right,key) in enumerate(zip([choices,selected,patients,stats],expected,keys)):
  if index in [1,2]:right=right[right.stage.ne('marker_only')]
  numerical.equal(left,right,key)
 tied=choices[choices.space.eq('noDR')&choices.method.eq('KMeans')]
 assert len(tied)==20 and tied.selected_k.eq(5).all()
 assert stats.n_patients_Lfine_ineligible.eq(1).all()
 # Perturb only one held-out fold. That fold's K choice must remain unchanged.
 target_fold=int(folds.fold.min());perturbed=data.copy()
 held=set(folds[folds.fold.eq(target_fold)].patient)
 mask=perturbed.patient.isin(held)&perturbed.k.eq(40)&perturbed.stage.eq('terminal090')&perturbed.lfine_metric_eligible
 perturbed.loc[mask,'lfine_macroF1']=.999
 alternative=numerical.select(perturbed,anchor,eligibility,folds,spec,protocol)[0]
 numerical.equal(choices[choices.fold.eq(target_fold)],alternative[alternative.fold.eq(target_fold)],keys[0])
 text=wrapper.compact_text(stats)
 assert text.count('\n| ')==22 and 'marker_only' not in text and 'L1' not in text
 OUT.mkdir(parents=True)
 nb=ROOT/'notebooks/dgscrna_results.ipynb';before=sha(nb)
 (OUT/'manifest.json').write_text(json.dumps(dict(status='synthetic_fixture_only',real_results=False)))
 wrapper.update_notebook(OUT,text,False)
 assert sha(nb)==before
 receipt=json.loads((OUT/'notebook_receipt.json').read_text());assert receipt['status']=='candidate' and receipt['remaining_cells_preserved'] and not receipt['PTC_changed']
 report=dict(status='passed_synthetic_only',real_campaign_aggregation_executed=False,notebook_unchanged=True,
  n_candidate_terminal_rows=len(data),n_anchor_terminal_rows=len(anchor),n_selection_rows=len(choices),n_contrasts=len(stats),
  original_three_stage_terminal_parity=True,independent_selector_parity=True,smaller_K_exact_ties=True,
  heldout_labels_cannot_change_own_K=True,truth_only_no_eligible_patient_retained=True,
  notebook_candidate_schema_and_preservation_passed=True,job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),
  sources={name:sha(CODE/name) for name in ['run.py','numerical.py','test_synthetic.py','SOURCE_PROTOCOL_MANIFEST.json']})
 (OUT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report,indent=2))
if __name__=='__main__':main()
