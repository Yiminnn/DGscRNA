from pathlib import Path
import json
import os
from ptc_common import BASE,RECOVERY,require_slurm,sha,write_json,utc
from label_rules import simplify,broad_lineage,strict_T,legacy_broad_T_sensitivity,T_subtype,MODULES,MAPPING_NOTES
require_slurm()
import pandas as pd
source=RECOVERY/'inventory_r/object_01_full_marker_symbol.content.json'
libs=json.loads(source.read_text())
rows=[]
for library,panels in libs.items():
 for name in panels:
  rows.append(dict(library=library,native=name,general=simplify(name),broad=broad_lineage(name),
       strict_T=strict_T(name),legacy_broad_T=legacy_broad_T_sensitivity(name),T_subtype=T_subtype(name)))
dest=BASE/'protocol';dest.mkdir(exist_ok=True)
pd.DataFrame(rows).to_csv(dest/'native_label_ontology.csv',index=False)
write_json(dest/'RNA_modules.json',MODULES)
write_json(dest/'label_rules_manifest.json',dict(status='frozen_before_evaluation',created_at=utc(),
   rule_source_sha256=sha(Path(__file__).with_name('label_rules.py')),
   libraries_source_sha256=sha(source),ontology_sha256=sha(dest/'native_label_ontology.csv'),
   modules_sha256=sha(dest/'RNA_modules.json'),notes=MAPPING_NOTES,job=os.environ['SLURM_JOB_ID'],
   annotation_predictions_read=False,TCR_outcomes_read=False))
(dest/'LABEL_RULES_FROZEN').write_text(sha(dest/'label_rules_manifest.json')+'\n')
print('Frozen name-only ontology for',len(rows),'library/panel entries')
