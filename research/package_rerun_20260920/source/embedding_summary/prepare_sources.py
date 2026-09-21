"""Copy frozen A1 Lfine selection protocol; preserve its numerical bodies."""
from pathlib import Path
import ast, hashlib, json
CODE=Path(__file__).resolve().parent
ROOT=CODE.parents[2]
OLD=ROOT/'handoff/reviewer_completion_20260920/embedding_lfine_v1'
CAMP=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
(CODE/'original').mkdir(exist_ok=True)
(CODE/'protocol').mkdir(exist_ok=True)
sources={}
for name in ['select_lfine.py','verify_selection.py','common.py']:
 source=OLD/name;(CODE/'original'/name).write_bytes(source.read_bytes());sources[str(source)]=sha(source)
protocols={'lfine.json':CAMP/'embedding_lfine_v1/protocol.json',
 'embedding.json':CAMP/'protocol/embedding.json',
 'selection.json':CAMP/'protocol/embedding_selection.json',
 'patient_folds.csv':CAMP/'protocol/embedding_patient_folds.csv'}
for name,source in protocols.items():
 (CODE/'protocol'/name).write_bytes(source.read_bytes());sources[str(source)]=sha(source)
common=(OLD/'common.py').read_text()
measures=next(ast.literal_eval(n.value) for n in ast.parse(common).body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='MEASURES' for t in n.targets))
source=(OLD/'select_lfine.py').read_text()
body=source[source.index('    selected_frames=[];choices=[];patients_all=[];paired=[];rng='):source.index('    out.mkdir(parents=True,exist_ok=True)')]
substitutions=[('c.MEASURES','MEASURES'),('n_samples*2*7*3*3','n_samples*2*7*3*2')]
for before,after in substitutions:
 assert before in body;body=body.replace(before,after)
header='"""Unchanged A1 patient selection/statistics; two terminal endpoints only."""\n'
header+='import numpy as np\nimport pandas as pd\nfrom scipy.stats import wilcoxon\nMEASURES='+repr(measures)+'\n\n'
text=header+'def select(allmetrics, anchors, eligibility, folds, spec, protocol):\n    foldmap=folds.drop_duplicates("patient").set_index("patient").fold.to_dict()\n'+body
text+='    return pd.DataFrame(choices), pd.concat(selected_frames), pd.concat(patients_all), stats\n\n'
source=(OLD/'verify_selection.py').read_text()
body2=source[source.index('    choices=[];selected_frames=[];patient_frames=[];comparisons=[];rng='):source.index('    equal(pd.DataFrame(choices),')]
substitutions2=[('c.MEASURES','MEASURES'),('n_samples*126','n_samples*84')]
for before,after in substitutions2:
 assert before in body2;body2=body2.replace(before,after)
text+='def independently_select(data, anchors, eligibility, folds, spec, protocol):\n    foldmap=folds.drop_duplicates("patient").set_index("patient").fold.to_dict()\n'+body2
text+='    return pd.DataFrame(choices), pd.concat(selected_frames), pd.concat(patient_frames), stats\n\n'
node=next(n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name=='equal')
text+=ast.get_source_segment(source,node)+'\n'
compile(text,str(CODE/'numerical.py'),'exec')
(CODE/'numerical.py').write_text(text)
manifest=dict(status='frozen_scientific_selection_body',source_hashes=sources,
 adapted_numerical_sha256=sha(CODE/'numerical.py'),measures=measures,
 substitutions={'selection':substitutions,'independent_verification':substitutions2},
 preserved='Training-patient selection, exact ties, bootstrap iteration order, Wilcoxon and Holm arithmetic unchanged',
 endpoint_difference='Only terminal070 and terminal090 retained; marker-only diagnostics omitted. Selection remains terminal090.',
 independent_inference=dict(path=str(CODE/'canonical_inference.py'),sha256=sha(CODE/'canonical_inference.py'),
  policy='Independently verify patient means/non-p statistics at original precision, then independently recompute Wilcoxon/Holm on lossless canonical producer vectors rebuilt from selected sample rows; no producer change or tolerance relaxation'),
 planned_counts=dict(representations=1694,candidate_threshold_rows=44044,anchor_threshold_rows=484,selection_rows=420,paired_contrasts=84))
(CODE/'SOURCE_PROTOCOL_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
