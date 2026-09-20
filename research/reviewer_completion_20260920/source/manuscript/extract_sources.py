#!/usr/bin/env python3
"""Read original document XML and record supplement containers; no data analysis."""
import hashlib
import json
import os
from pathlib import Path
import re
import zipfile
import xml.etree.ElementTree as ET

assert os.environ.get('SLURM_JOB_ID'), 'Use SLURM for document extraction'
ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT = ROOT / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/manuscript'
OUT.mkdir(parents=True, exist_ok=True)
NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
documents = {
    'original_v16': ROOT/'paper/submission_v16/DG_scRNA_04232026_V16_cell_report.docx',
    'existing_reply': ROOT/'paper/revision/Review_Reply.docx',
    'original_delivery_table': ROOT/'results/hvg_ptc_20260916_v1/geo_submission_v1/PTC_GEO_submission_20260916/05_provenance/original_delivery_sample_table.docx',
}

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()

records = []
for name, source in documents.items():
    with zipfile.ZipFile(source) as archive:
        tree = ET.fromstring(archive.read('word/document.xml'))
    paragraphs = []
    for index, paragraph in enumerate(tree.findall('.//w:p', NS), 1):
        text = ''.join(node.text or '' for node in paragraph.findall('.//w:t', NS))
        if text.strip():
            paragraphs.append(f'P{index:04d}\t{text}')
    target = OUT/f'{name}.paragraphs.txt'
    target.write_text('\n'.join(paragraphs)+'\n')
    records.append({'source': str(source.relative_to(ROOT)), 'sha256': sha(source),
                    'extraction': str(target.relative_to(ROOT)), 'nonempty_paragraphs': len(paragraphs),
                    'note': 'w:t text; table paragraphs included; paragraph IDs are XML positions, not page numbers'})

supplements = []
for path in sorted((ROOT/'paper/submission_v16').glob('Supplementary*')):
    item = {'path': str(path.relative_to(ROOT)), 'bytes': path.stat().st_size, 'sha256': sha(path)}
    if path.suffix == '.xlsx':
        with zipfile.ZipFile(path) as archive:
            bad = archive.testzip()
            item['zip_crc_integrity'] = 'passed' if bad is None else bad
            n = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
            book = ET.fromstring(archive.read('xl/workbook.xml'))
            item['sheet_names'] = [e.attrib['name'] for e in book.findall('.//s:sheet', n)]
            item['sheets_xml_dimensions'] = {}
            for member in archive.namelist():
                if re.fullmatch(r'xl/worksheets/sheet\d+\.xml', member):
                    with archive.open(member) as stream:
                        first = stream.read(4096).decode('utf-8')
                    dim = re.search(r'<dimension ref="([^"]+)"', first)
                    item['sheets_xml_dimensions'][member] = dim.group(1) if dim else 'not in XML prefix'
        item['remaining'] = 'Full spreadsheet reader, endpoint/content validation and submission cross-reference check remain; XML dimensions are metadata, not recomputed cell counts.'
    else:
        item['remaining'] = 'Final figure legibility/caption/numbering and reader rendering need verification.'
    supplements.append(item)

(OUT/'source_inventory.json').write_text(json.dumps({'slurm_job_id': os.environ['SLURM_JOB_ID'],
    'slurm_step_id': os.environ.get('SLURM_STEP_ID'), 'documents': records,
    'original_supplements': supplements}, ensure_ascii=False, indent=2)+'\n')
print(json.dumps({'documents': len(records), 'supplement_files': len(supplements)}))
