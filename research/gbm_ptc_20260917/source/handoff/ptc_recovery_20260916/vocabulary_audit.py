"""Report marker-library vocabulary capacity without changing any annotation."""
import json
import os
from pathlib import Path
from ptc_common import BASE,require_slurm,sha,utc,write_json
require_slurm()
import pandas as pd
ontology=pd.read_csv(BASE/'protocol/native_label_ontology.csv')
ref=pd.read_csv(BASE/'evaluation_reference/reference_cells.csv.gz',usecols=['sample','patient','group','S2_terminal_broad','S3_terminal_broad'])
rows=[];libraries=[]
for library,panels in ontology.groupby('library'):
    supported=set(panels.broad)|{'Unknown'}
    libraries.append(dict(library=library,n_native_panels=len(panels),n_distinct_broad_lineages=panels.broad.nunique(),
       contains_strict_T='T' in supported,broad_lineages=json.dumps(sorted(set(panels.broad)))))
    for sample,r in ref.groupby('sample'):
      for source in ['S2','S3']:
        y=r[source+'_terminal_broad'];present=set(y)
        rows.append(dict(library=library,sample=sample,patient=r.patient.iloc[0],group=r.group.iloc[0],reference=source,
           n_cells=len(y),n_reference_broad_classes=len(present),n_supported_reference_classes=len(present&supported),
           theoretical_macro_F1_present_ceiling=len(present&supported)/len(present),
           reference_cell_fraction_in_supported_vocabulary=y.isin(supported).mean(),
           unsupported_reference_classes=json.dumps(sorted(present-supported)),
           supports_Unknown_abstention=True))
out=BASE/'vocabulary_audit';out.mkdir(exist_ok=True)
pd.DataFrame(libraries).to_csv(out/'original_library_vocabulary.csv',index=False)
pd.DataFrame(rows).to_csv(out/'reference_vocabulary_capacity_by_sample.csv',index=False)
write_json(out/'manifest.json',dict(status='complete',n_libraries=17,n_reference_sample_rows=len(rows),
    source_sha256=sha(Path(__file__)),ontology_sha256=sha(BASE/'protocol/native_label_ontology.csv'),
    interpretation='Name-based representational capacity only, not achieved accuracy. Assumes perfect calls within supported labels and Unknown abstention elsewhere. Actual marker detection, seed availability and refinement can lower performance. S2/S3 are historical outputs, not independent truth.',
    prediction_labels_changed=False,outputs={p.name:sha(p) for p in out.glob('*.csv')},
    completed_at=utc(),job=os.environ['SLURM_JOB_ID']))
(out/'COMPLETE').write_text(sha(out/'manifest.json')+'\n')
print('Original 17-library vocabulary and reference-class capacity audited',flush=True)
