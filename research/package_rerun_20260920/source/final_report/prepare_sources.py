"""Freeze old B rendering and C selected-patient reduction sources only."""
from pathlib import Path
import ast,hashlib,json
CODE=Path(__file__).resolve().parent;ROOT=CODE.parents[2]
OLD=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
sources={}
source=ROOT/'handoff/reviewer_completion_20260920/controls_lfine_presentation_v1/build.py'
text=source.read_text();sources[str(source)]=sha(source)
body=text[text.index("samples = ['TKU4163'"):text.index('counts = pd.DataFrame(descriptions)')]
header='"""Exact original B plotting/checkpoint body with caller-supplied fresh rows."""\nimport json\nimport numpy as np\nimport pandas as pd\nimport matplotlib\nmatplotlib.use("Agg")\nimport matplotlib.pyplot as plt\nfrom matplotlib.lines import Line2D\n'
code=header+'\ndef render(frame, OUT):\n    assert len(frame)==408 and frame.row_key.is_unique and frame.terminal_valid.all()\n    OUT.mkdir(parents=True)\n'+''.join('    '+line+'\n' for line in body.splitlines())+'    return pd.DataFrame(descriptions), checkpoints\n'
compile(code,str(CODE/'render_B.py'),'exec');(CODE/'render_B.py').write_text(code)
canonical=ROOT/'handoff/reviewer_completion_20260920/comparison_lfine_v1/selection_v3/canonical.py'
text=canonical.read_text();sources[str(canonical)]=sha(canonical)
names={'ordered_patient_means','ordered_training_rank'};functions=[];constants=[]
for node in ast.parse(text).body:
 if isinstance(node,ast.FunctionDef) and node.name in names:functions.append(ast.get_source_segment(text,node))
 if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id in {'SPEC','MEASURES'} for t in node.targets):constants.append(ast.get_source_segment(text,node))
(CODE/'canonical_C.py').write_text('"""Copied C ordered reductions; no original imports or execution."""\n'+'\n\n'.join(constants+functions)+'\n')
metadata=[OLD/'controls_lfine_v1/full/validation.json',OLD/'controls_lfine_v1/independent_verification/validation.json',
 OLD/'controls_lfine_presentation_v1/manifest.json',OLD/'controls_lfine_presentation_v1/validation.json',
 OLD/'no_clustering_lfine_v1/summary/manifest.json',OLD/'no_clustering_lfine_v1/full_validation/validation.json',
 OLD/'no_clustering_lfine_v1/protocol.json',OLD/'no_clustering/patient_folds.csv',
 OLD/'comparison_lfine_v1/summary/manifest.json',OLD/'comparison_lfine_v1/selection_validation/validation.json',
 ROOT/'handoff/reviewer_completion_20260920/comparison_lfine_v1/selection_v3/selection_freeze.json',
 ROOT/'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/protocol/patient_folds.csv']
for path in metadata:sources[str(path)]=sha(path)
for path in [ROOT/'handoff/reviewer_completion_20260920/no_clustering_lfine_v1/aggregate.py',
 ROOT/'handoff/reviewer_completion_20260920/comparison_lfine_v1/selection_v3/select_patients.py',
 ROOT/'handoff/reviewer_completion_20260920/compact_B_A2_handoff_v1/prepare.py']:
 sources[str(path)]=sha(path)
manifest=dict(status='source_preparation_only',source_hashes=sources,
 derived_files={name:sha(CODE/name) for name in ['render_B.py','canonical_C.py']},
 scientific_derivation='B plotting/checkpoint body copied unchanged from samples declaration through REPORT write; C ordered means/rank functions and constants copied unchanged by AST.',
 retained_inference='A2 and non-DG C confidence intervals/p values retained only after all relevant fresh terminal metric vectors and fixed selection mappings reproduce their frozen inputs; no new inference claimed.')
(CODE/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
