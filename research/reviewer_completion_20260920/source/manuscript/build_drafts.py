#!/usr/bin/env python3
"""Assemble a source-linked writing draft; no experiments or computed metrics."""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
CODE = ROOT/'handoff/reviewer_completion_20260920/manuscript'
OUT = ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/manuscript'
PLAN = ROOT/'handoff/reviewer_completion_plan_20260920'
OUT.mkdir(parents=True, exist_ok=True)
rows = list(csv.DictReader((PLAN/'REVIEWER_COVERAGE.tsv').open(), delimiter='\t'))
content = json.loads((CODE/'response_content.json').read_text())
coverage = {row['id']: row for row in json.loads((PLAN/'plan_coverage.json').read_text())['rows']}
assert len(rows) == len(content) == 35 and {row['id'] for row in rows} == set(content)
inventory = json.loads((OUT/'source_inventory.json').read_text())
parts = ['''# Response to reviewers — working revision draft

Date: 2026-09-20. **Internal draft; not a submitted or completed response.**

This document covers all 35 requirements in the approved completion plan. It
preserves the original manuscript and the older response document. The response
paragraphs below are proposed wording. A statement about future revision is not
evidence that the final manuscript has been edited. Bracketed evidence slots must
be replaced with validated artifacts and final page/figure/table locations before
submission. Running jobs, completed jobs and scientifically accepted results are
distinct states. No new numerical result is claimed by this drafting task.

The reviewer source is `paper/comments.md`, dated 2026-07-09. Reviewer IDs follow
`handoff/reviewer_completion_plan_20260920/REVIEWER_COVERAGE.tsv`; the older July
reply uses a different subdivision and must not supply the current ID mapping.
The source manuscript is
`paper/submission_v16/DG_scRNA_04232026_V16_cell_report.docx`; source extraction,
hashes and supplement-container evidence are in `source_inventory.json`.

Fixed historical policy: preserve the original SignacX None / zero-predicted-T
result and original S3 labels. Preserve scType's actual result. Later nonzero
SignacX runs remain separate and are neither hidden nor forced to zero.
''']
pending_rows = []
for row in rows:
    identifier = row['id']
    title, paragraph, pending = content[identifier]
    refs = coverage[identifier].get('evidence_paths', [])
    parts += [f'## {identifier} — {title}\n',
              f'Reviewer source: `{row["source"]}`. Work packages: {row["work_packages"]}.\n',
              '**Proposed response**\n\n' + paragraph + '\n',
              '**Pending acceptance slot**\n\n[' + pending + ']\n',
              '**Evidence register** — these are existing source/analysis locations, not a declaration that all current acceptance criteria are met.\n']
    for ref in refs:
        parts.append(f'- `{ref}`')
    parts.append('')
    pending_rows.append({'id': identifier, 'response_title': title,
        'work_packages': row['work_packages'], 'source': row['source'],
        'draft_status': 'drafted_pending_evidence_and_manuscript_integration',
        'acceptance_slot': pending, 'final_manuscript_location': 'PENDING',
        'final_figure_table': 'PENDING'})
parts.append('''## Finalization gate

Before any response is described as complete, replace every acceptance slot with
verified evidence or an explicit source-supported limitation and the corresponding
manuscript change. Check all35 IDs against the approved acceptance table, including
external-source requirements for wet-lab facts, raw reads, actual GEO access and
permanent archival. Produce and inspect the final submission PDF and supplements.
Do not convert partial computational coverage into a claim that all reviewer
requirements have been closed.
''')
(OUT/'RESPONSE_DRAFT_EN.md').write_text('\n'.join(parts))
with (OUT/'RESPONSE_ACCEPTANCE_SLOTS.tsv').open('w') as stream:
    writer = csv.DictWriter(stream, delimiter='\t', fieldnames=list(pending_rows[0]))
    writer.writeheader();writer.writerows(pending_rows)

supp = ['''# Original supplements: inventory and remaining verification

Inventory date: 2026-09-20. Source files are preserved in `paper/submission_v16/`.
Document extraction and container checks ran in SLURM7204040.43. XLSX CRC and
workbook-XML checks establish container readability, not correct biological
content or final spreadsheet-reader rendering. No matrix loading or new PTC
computation occurred in this documentation task.

The original V16 supplement descriptions are in extracted paragraphs P0385-P0389.
New Supplement A/B must supplement these requested original files, not silently
replace their identity or endpoints.
''']
for item in inventory['original_supplements']:
    supp.append(f'## {Path(item["path"]).name}\n')
    supp.append(f'Original path: `{item["path"]}`\n\nSHA-256: `{item["sha256"]}`\n\nBytes: {item["bytes"]}.\n')
    if 'sheet_names' in item:
        supp.append('Container CRC: ' + item['zip_crc_integrity'] + '. Workbook sheets: ' + '; '.join('`'+s+'`' for s in item['sheet_names']) + '.\n')
        supp.append('Declared XML sheet dimensions: ' + '; '.join(f'`{k}: {v}`' for k,v in item['sheets_xml_dimensions'].items()) + '. These are metadata, not newly computed cell counts.\n')
    supp.append('Remaining: ' + item['remaining'] + '\n')
supp.append('''## Content-specific acceptance

| Original item | Required revision checks |
|---|---|
| Figure S1 | Render the original PDF/PNG and revised submission figure; verify legible labels, method names, preprocessing/input provenance and captions. A separated UMAP picture does not establish superior accuracy or preserved original-space density. Bind any empirical claim to actual partition/terminal metrics. |
| Table S1 | Check gene lists, tissue/context and assay provenance; distinguish AllTissues seven-tissue union from all-human/all-tissue catalogues. The first workbook sheet is named `Supplementary Table 2` despite the S1 filename; correct only in a versioned revision copy and document the change. Map this table to the exact marker configurations used in current runs. |
| Table S2 | Preserve the original final annotations; verify sample/barcode identity, headers and endpoint. Do not rename TCR support as independent T-cell-subtype validation. Reconcile its historical labels independently of S3. |
| Table S3 | Preserve native predictions, original flags and historical SignacX zero-T row. Separate the paper Table2 metric reconstruction, literal S3 evaluation and new strict-T comparisons. Verify all-cell denominators, ontology mappings and historical Accuracy disclosure. |
| Table S4 | Trace each composition value to its S2/S3 endpoint, full denominator and sample/patient identity; retain paired patient structure and contextual limits. |

Final delivery needs a supplement index, direct readable previews where useful,
machine-readable tables, consistent numbering and tested local/packaged links.
The current inventory does not claim this final submission-level acceptance.
''')
(OUT/'SUPPLEMENT_INVENTORY.md').write_text('\n'.join(supp))

for filename in ['METHODS_REPLACEMENT.md','SOURCE_GAPS.md','response_content.json','extract_sources.py','build_drafts.py']:
    shutil.copy2(CODE/filename, OUT/filename)
for filename in ['RESPONSE_DRAFT_EN.md','RESPONSE_ACCEPTANCE_SLOTS.tsv','SUPPLEMENT_INVENTORY.md']:
    shutil.copy2(OUT/filename, CODE/filename)

readme = '''# Reviewer G/H drafting package

This is a source-grounded documentation milestone, not completed experiments,
not a rewritten final submission, and not acceptance of entire work packages G/H.

- `RESPONSE_DRAFT_EN.md`: 35-point English draft with explicit pending slots.
- `RESPONSE_ACCEPTANCE_SLOTS.tsv`: all35 IDs and unresolved final locations.
- `METHODS_REPLACEMENT.md`: concrete original-R-aligned replacement wording.
- `SOURCE_GAPS.md`: wet-lab/software/chemistry/count-flow ledger and source gaps.
- `SUPPLEMENT_INVENTORY.md`: original FigureS1 and S1–S4 paths, hashes and checks.
- `source_inventory.json`: immutable original-document hashes and container check.
- `original_v16.paragraphs.txt`, `existing_reply.paragraphs.txt`,
  `original_delivery_table.paragraphs.txt`: source extracts with XML paragraph IDs.

Original files were not modified. Document extraction and CRC inspection used
SLURM7204040.43. No PTC computation, external publication, connector or account
access occurred. Bibliography verification against external primary sources,
fresh-environment execution, final manuscript/PDF editing and final supplement
rendering remain separate tasks.
'''
(OUT/'README.md').write_text(readme)
(CODE/'README.md').write_text(readme)

status = {'stage':'G_H_DOCS','status':'draft_ready_for_review',
    'updated_at':datetime.now(timezone.utc).isoformat(),
    'summary':'35-point English draft, original-R Methods corrections, source-gap ledger and original-supplement inventory prepared; final manuscript and G/H acceptance remain pending.',
    'completed_documentation_items':['35 source-linked reviewer response draft entries','Methods replacement wording','wet-lab and cell-flow source-gap ledger','original Figure S1/S1-S4 container inventory'],
    'not_completed':['ongoing experiments and statistical acceptance','final manuscript/PDF revision','missing original wet-lab/run provenance','cell-level attrition reconstruction','final supplement reader/content validation','authorized GEO/DOI publication','clean-environment reproduction'],
    'jobs':[{'job_id':'7204040.43','purpose':'document XML extraction and original supplement CRC inventory','state':'completed'}],
    'evidence':[str((OUT/name).relative_to(ROOT)) for name in ['RESPONSE_DRAFT_EN.md','METHODS_REPLACEMENT.md','SOURCE_GAPS.md','SUPPLEMENT_INVENTORY.md','source_inventory.json']],
    'requirements':[row['id'] for row in rows],
    'validation':{'response_ids_match_approved_plan':True,'response_count':35,'new_results_claimed':False,'original_manuscript_modified':False,'full_work_packages_complete':False}}
(OUT/'status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2)+'\n')
(CODE/'status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2)+'\n')

manifest = {'status':'drafted_not_submission_complete',
    'files':[{'path':str(path.relative_to(OUT)),'bytes':path.stat().st_size,
              'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
             for path in sorted(OUT.iterdir()) if path.is_file() and path.name!='manifest.json']}
(OUT/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':'draft_ready_for_review','response_count':len(content),'output':str(OUT)}))
