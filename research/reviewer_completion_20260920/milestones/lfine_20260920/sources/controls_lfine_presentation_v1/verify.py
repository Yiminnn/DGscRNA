"""Independent B presentation-data audit; does not rebuild or alter artifacts."""
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
import hashlib,json,os,statistics
import numpy as np
import pandas as pd
from PIL import Image

assert os.environ.get('SLURM_JOB_ID')
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
CAMP=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
OUT=CAMP/'controls_lfine_presentation_v1';SRC=CAMP/'controls_lfine_v1'
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(Path(p).read_text())
def near(a,b):assert np.isclose(a,b,rtol=0,atol=1e-12,equal_nan=True),(a,b)

manifest=read(OUT/'manifest.json');original_hash=sha(OUT/'manifest.json')
assert (OUT/'COMPLETE').read_text().strip()==original_hash
assert manifest['source_sha256']==sha(ROOT/'handoff/reviewer_completion_20260920/controls_lfine_presentation_v1/build.py')
for path,digest in manifest['inputs'].items():assert sha(ROOT/path)==digest,path
for name,digest in manifest['files'].items():assert sha(OUT/name)==digest,name
proof=read(SRC/'independent_verification/validation.json');full=read(SRC/'full/validation.json')
assert proof['status']=='passed' and full['status']=='passed' and proof['n_rows']==408
assert proof['full_validation_sha256']==sha(SRC/'full/validation.json')
for name,digest in full['files'].items():assert sha(SRC/'full'/name)==digest
source=pd.read_csv(SRC/'full/terminal_metrics.csv.gz')
assert len(source)==408 and source.row_key.is_unique
assert source.terminal_valid.dtype==bool and source.terminal_valid.all()
neighbors=source[(source.task=='neighbors')&(source.stage=='terminal090')]
assert len(neighbors)==132
points=pd.read_csv(OUT/'neighbor_plot_points.csv.gz')
description=pd.read_csv(OUT/'neighbor_descriptive_summary.csv')
assert len(points)==144 and points.row_key.nunique()==132 and len(description)==48
assert set(points.row_key)==set(neighbors.row_key)
assert not points.duplicated(['panel','row_key']).any()
expected=neighbors.set_index('row_key')
pd.testing.assert_frame_equal(points[source.columns].reset_index(drop=True),
    expected.loc[points.row_key].reset_index()[source.columns],check_exact=False,rtol=0,atol=1e-12)
frequency=Counter(points.row_key)
assert Counter(frequency.values())=={1:120,2:12}
for rowkey,n in frequency.items():
    row=expected.loc[rowkey]
    shared=row.route=='UMAP2_SNN' and row.snn_k==20 and row.umap_neighbors==30
    assert (n==2)==shared
    if shared:assert set(points.loc[points.row_key==rowkey,'parameter'])=={'snn_k','umap_neighbors'}
families=points.groupby(['sample','budget','route','library','parameter'])
assert len(families)==48 and families.size().eq(3).all()
relations=Counter()
for (sample,budget,route,library,parameter),g in families:
    assert len(set(g.panel))==1
    setting=[10,20,40] if parameter=='snn_k' else [15,30,60]
    default=20 if parameter=='snn_k' else 30
    other='umap_neighbors' if parameter=='snn_k' else 'snn_k'
    assert sorted(g.x)==setting and np.array_equal(g.x,g[parameter])
    assert g[other].eq(30 if other=='umap_neighbors' else 20).all()
    assert route in (['PCA30_SNN','UMAP2_SNN'] if parameter=='snn_k' else ['UMAP2_SNN','UMAP2_HDBSCAN_R'])
    assert g.panel.eq(f'{budget}/{route}/{parameter}').all()
    d=description[(description['sample']==sample)&(description.budget==budget)&(description.route==route)&(description.library==library)&(description.parameter==parameter)]
    assert len(d)==1;d=d.iloc[0];baseline=g[g.x==default].iloc[0]
    assert d.default==default
    near(d.default_lfine_macroF1,baseline.lfine_macroF1);near(d.default_coverage,baseline.coverage)
    near(d.tested_min_lfine_macroF1,min(g.lfine_macroF1));near(d.tested_max_lfine_macroF1,max(g.lfine_macroF1))
    maximum=max(g.lfine_macroF1)
    rel='lower_than_an_alternative' if maximum>baseline.lfine_macroF1+1e-12 else ('tied_maximum' if sum(abs(v-maximum)<=1e-12 for v in g.lfine_macroF1)>1 else 'unique_maximum')
    assert d.default_relation==rel;relations[library,rel]+=1
    diagnostic=json.loads(d['values']);assert len(diagnostic)==3
    for saved in diagnostic:
        r=g[g.row_key==saved['row_key']];assert len(r)==1;r=r.iloc[0]
        assert saved['setting']==r.x;near(saved['lfine_macroF1'],r.lfine_macroF1);near(saved['coverage'],r.coverage)
assert relations==Counter({(r['library'],r['relation']):r['n'] for r in manifest['descriptive_default_relations']})
learning=source[source.task=='learning'];checkpoint=pd.read_csv(OUT/'checkpoint_lfine_all_thresholds.csv')
assert len(learning)==144 and len(checkpoint)==48
assert set(learning.stage)==set(checkpoint.stage)=={'terminal090','terminal070'}
assert set(learning.epochs)=={5,10,20,30} and set(learning.model_seed)=={0,1,42}
assert set(learning.library)=={'CM2_glioma_other'}
assert not checkpoint.duplicated(['sample','budget','stage','epochs']).any()
groups=learning.groupby(['sample','budget','stage','epochs'])
for keys,g in groups:
    assert len(g)==3 and set(g.model_seed)=={0,1,42}
    sample,budget,stage,epochs=keys
    row=checkpoint[(checkpoint['sample']==sample)&(checkpoint.budget==budget)&(checkpoint.stage==stage)&(checkpoint.epochs==epochs)]
    assert len(row)==1;row=row.iloc[0]
    # Explicit sample-SD denominator n-1, independently of pandas aggregation.
    for source_name,target in [('lfine_macroF1','lfine_macroF1'),('coverage','coverage')]:
        values=g[source_name].tolist();mean=sum(values)/3
        sd=(sum((v-mean)**2 for v in values)/2)**.5
        near(row[target+'_mean'],mean);near(row[target+'_seed_SD'],sd)
    assert row.n_seeds==3 and row.known_training_classes==min(g.n_training_classes)
with Image.open(OUT/'neighbor_lfine.png') as png:
    png.verify()
with Image.open(OUT/'neighbor_lfine.png') as png:dimensions=png.size
assert (OUT/'neighbor_lfine.pdf').read_bytes().startswith(b'%PDF-')
report=read(SRC/'full/validation.json')
assert sha(OUT/'manifest.json')==original_hash
assert manifest['n_checkpoint_aggregate_rows']==48 and manifest['n_descriptive_sweep_families']==48
result=dict(status='passed',scope='Independent B Lfine presentation data/source audit; root visual review remains separate',
    presentation_manifest_sha256=original_hash,validation_source_sha256=sha(Path(__file__)),
    source_validation_sha256=sha(SRC/'independent_verification/validation.json'),
    n_source_rows=408,n_unique_neighbor_terminal090_rows=132,n_plotted_neighbor_points=144,
    n_shared_UMAP_SNN_default_points=12,n_sweep_families=48,n_checkpoint_seed_rows=144,n_checkpoint_aggregate_rows=48,
    all_neighbor_rows_and_all_source_fields_equal=True,all48_default_relation_descriptions_recomputed=True,
    all_checkpoint_means_sample_SDs_and_coverages_recomputed=True,thresholds=[.7,.9],epochs=[5,10,20,30],model_seeds=[0,1,42],
    learning_coverage_equals_legacy_called_coverage=bool(np.allclose(learning.coverage,learning.legacy_called_coverage,rtol=0,atol=1e-12)),
    source_hashes_verified=True,PNG_dimensions=list(dimensions),PDF_header_valid=True,
    graphical_scope='Saved plot-point data independently verified; frozen builder asserts each line x/y equals those inputs; no independent re-render or claim of visual review',
    no_fitting=True,original_artifacts_unchanged=True,no_canonical_notebook_writes=True,
    job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),completed_at=datetime.now(timezone.utc).isoformat())
dest=OUT/'validation.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2),flush=True)
