"""Independent read-only notebook preservation, chapter and saved-figure audit."""
from pathlib import Path
from datetime import datetime, timezone
import base64
import csv
import hashlib
import json
import os
import re
import nbformat

assert os.environ.get('SLURM_JOB_ID')
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
CAMP=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
OUT=CAMP/'notebook'
TAG='reviewer_completion_20260920'
TARGET=ROOT/'notebooks/dgscrna_results.ipynb'
BACKUP=OUT/'before_9c38f579348b1052.ipynb'


def sha(path):
    with path.open('rb') as handle:return hashlib.file_digest(handle,'sha256').hexdigest()


before_hash=sha(BACKUP);after_hash=sha(TARGET)
manifest=json.loads((OUT/'manifest.json').read_text())
assert before_hash==manifest['input_notebook_sha256']=='9c38f579348b1052523cfdae839847389daf5e7b5f7b0eb32c12e66183214afb'
assert after_hash==manifest['notebook_sha256']
before=nbformat.read(BACKUP,as_version=4);after=nbformat.read(TARGET,as_version=4)
nbformat.validate(after)
old=[cell for cell in after.cells if TAG not in cell.metadata.get('tags',[])]
new=[cell for cell in after.cells if TAG in cell.metadata.get('tags',[])]
assert len(before.cells)==len(old)==1067 and len(new)==13 and len(after.cells)==1080
assert old==before.cells, 'An original cell/source/output/metadata or original ordering changed'
assert after.metadata==before.metadata, 'Notebook-level metadata changed'
assert (after.nbformat,after.nbformat_minor)==(before.nbformat,before.nbformat_minor)
ptc_indices=[i for i,cell in enumerate(before.cells) if re.search(r'\bPTC\b|thyroid|\bTU-[12]\b|\bMT-[12]\b',cell.source,re.I)]
assert all(old[i]==before.cells[i] for i in ptc_indices)
expected_counts={'1':1,'4':5,'5':1,'9':6};actual_counts={key:0 for key in expected_counts}
last_original={chapter:max(i for i,cell in enumerate(before.cells) if cell.metadata.get('v5_option_a_chapter')==chapter) for chapter in expected_counts}
prior_original=-1
for cell in after.cells:
    if TAG not in cell.metadata.get('tags',[]):
        prior_original+=1
    else:
        chapter=cell.metadata['v5_option_a_chapter']
        assert chapter in expected_counts
        assert prior_original==last_original[chapter], 'Addition is outside its original chapter boundary'
        assert cell.metadata['v5_option_a_role']=='reviewer_addition'
        actual_counts[chapter]+=1
assert actual_counts==expected_counts
errors=[output for cell in new for output in cell.get('outputs',[]) if output.get('output_type')=='error']
assert not errors
expected_images={
    'evidence/NL022_selected_marker_dotplot.png':'4',
    'evidence/NL022_selected_marker_violin.png':'4',
    'evidence/NL022_markers_initial_terminal_same_coordinates.png':'4',
    'embedding/TKU4163/hvg2000/PCA2/figures/KMeans_K10.png':'9',
    'embedding/TKU4163/hvg2000/PCA2/figures/HDBSCAN_R.png':'9',
}
image_records=[]
for cell in new:
    encoded=[output.get('data',{}).get('image/png') for output in cell.get('outputs',[]) if 'image/png' in output.get('data',{})]
    if not encoded:continue
    match=re.fullmatch(r"review_image\(REVIEW_CAMP/'([^']+)'\)",cell.source.strip())
    assert match and len(encoded)==1
    relative=match.group(1)
    assert relative in expected_images and cell.metadata['v5_option_a_chapter']==expected_images[relative]
    encoded=encoded[0]
    image_bytes=base64.b64decode(''.join(encoded) if isinstance(encoded,list) else encoded)
    digest=hashlib.sha256(image_bytes).hexdigest()
    path=CAMP/relative
    assert digest==sha(path)
    source_manifest=json.loads((path.parent/'manifest.json').read_text())
    assert digest==source_manifest['files'][path.name]
    image_records.append(dict(source=str(path),source_sha256=digest,embedded_PNG_exact=True,
                              chapter=cell.metadata['v5_option_a_chapter'],source_manifest_sha256=sha(path.parent/'manifest.json')))
assert len(image_records)==5 and {record['source'] for record in image_records}=={str(CAMP/path) for path in expected_images}
index=OUT/'all_completed_A1_figures.csv'
with index.open() as handle:records=list(csv.DictReader(handle))
assert len(records)==manifest['all_A1_figure_index_rows']==52
assert len({row['figure_relative_to_notebook'] for row in records})==52
for row in records:
    path=(TARGET.parent/row['figure_relative_to_notebook']).resolve()
    assert path.is_relative_to(ROOT) and sha(path)==row['sha256']
    fm=json.loads((path.parent/'manifest.json').read_text())
    assert fm['files'][path.name]==row['sha256']
joined='\n'.join(cell.source for cell in new if cell.cell_type=='markdown')
for required in ['Cohort rankings await complete results and verification',
                 'consistency rather than independent biological validation or optimality',
                 'corrected121-sample run is in progress',
                 'historical SignacX None / zero predicted T-cell entry is retained']:
    assert required in joined
assert sha(TARGET)==after_hash and sha(BACKUP)==before_hash, 'Notebook changed during read-only verification'
report=dict(status='passed',scope='Saved notebook snapshot only; full GBM experiments remain in progress',
            notebook=str(TARGET),notebook_sha256=after_hash,backup_sha256=before_hash,
            original_cells=1067,original_cells_source_output_metadata_and_order_exact=True,
            notebook_metadata_and_format_exact=True,ptc_cells_exact=True,ptc_source_match_cells=len(ptc_indices),
            inserted_cells=13,total_cells=1080,insertions_by_existing_chapter=actual_counts,new_output_errors=0,
            decoded_saved_PNGs=5,images=image_records,all_completed_A1_index_rows=52,
            all_indexed_figure_source_hashes_match=True,scientific_claims_mark_pending_and_descriptive=True,
            no_notebook_execution_or_model_fitting=True,validation_source_sha256=sha(Path(__file__)),
            job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),completed_at=datetime.now(timezone.utc).isoformat())
(OUT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
