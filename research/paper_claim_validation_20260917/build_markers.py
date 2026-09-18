"""Freeze previously curated and anatomy-selected marker libraries before scoring."""
import json
import os
import re
import shutil
from pathlib import Path
from common import ROOT, OUT, OLD, L1, require_slurm, sha, write_json, complete, utc

def native_term(panel):
    if panel.startswith('CM2+'):
        return '+'.join(panel.split('+')[4:])
    return panel.split('+')[-1]

def semantic(term):
    """Fixed broad taxonomy, never use a sample's labels to map predictions."""
    s=term.lower().replace('‐','-').replace('_',' ').strip()
    s=re.sub(r'\s+',' ',s)
    if s in ['unknown','undecided','no annotation','']:return 'Unknown'
    if 'non-neuron' in s or 'non neuron' in s:return 'UNMAPPABLE'
    if any(t in s for t in ['cancer cell','cancer stem','glioma stem','glioblastoma cell','malignant','neoplastic']):return 'Malignant'
    if any(t in s for t in ['oligodendrocyte precursor','oligodendrocyte progenitor']) or s=='opc':return 'OPC'
    if 'oligodendrocyte' in s:return 'Oligodendrocyte'
    if 'astrocyte' in s:return 'Astrocyte'
    if any(t in s for t in ['inhibitory neuron','interneuron','gabaergic','somatostatin interneuron']):return 'Inhibitory neuron'
    if any(t in s for t in ['excitatory neuron','glutamatergic neuron','pyramidal neuron']):return 'Excitatory neuron'
    if 'neuron' in s:return 'AMBIGUOUS_NEURON'
    if any(t in s for t in ['macrophage','microglia','monocyte','myeloid','dendritic','neutrophil']):return 'TAM'
    if (any(t in s for t in ['lymphocyte','natural killer','plasma cell','plasmablast','plasmocyte'])
        or re.search(r'\b[bt][ -]cell',s) or re.search(r'\bt\(',s)):
        return 'Lymphocyte'
    if 'endothel' in s:return 'Endothel'
    if any(t in s for t in ['pericyte','mural cell','vascular smooth muscle']):return 'Pericyte'
    return 'UNMAPPABLE'

def run():
    require_slurm()
    import pandas as pd
    d=OUT/'markers';d.mkdir(parents=True,exist_ok=True)
    if (d/'MARKERS_COMPLETE').exists():
        assert (d/'MARKERS_COMPLETE').read_text().strip()==sha(d/'manifest.json');return
    oldmap=pd.read_csv(ROOT/'handoff/markers_v3/mapping_L1_v3.csv',keep_default_na=False)
    olddict={(r.marker_set,r.panel):r.gold_class for r in oldmap.itertuples()}
    prior_terms={}
    for r in oldmap.itertuples():
        if r.panel.startswith('CM2_'):
            term='+'.join(r.panel.split('+')[2:]).strip().lower()
            prior_terms.setdefault(term,set()).add(r.gold_class)
    libraries={};sources=[];maps=[]
    # First library is the historical fixed GSE274546 comparator, not newly selected.
    names=['CM2_glioma_other','CM2_glioblastoma','CM2_brain_normal','CARE_TME',
           'BrainAtlas112','Liu_devbrain','UNION_CellMarker2_brain','UNION_all']
    for name in names:
        src=ROOT/'handoff/markers_v3'/f'{name}.csv'
        t=pd.read_csv(src,index_col=0,keep_default_na=False,dtype=str)
        libraries[name]={p:list(dict.fromkeys(g.strip() for g in t[p] if g.strip())) for p in t}
        assert all(libraries[name].values())
        sources.append(dict(library=name,source=str(src),sha256=sha(src),
            type='previously curated marker library',
            reference_label_construction_overlap=name in ['CARE_TME','BrainAtlas112','UNION_all'],
            evidence_audit='CARE_TME and BrainAtlas112 were used in author reference-label construction; UNION_all includes them. Their results are concordance, not independent validation. Other panels also require source-study/species/assay audit.'))
        for panel in t:
            value=olddict.get((name,panel),'UNMAPPABLE')
            # Preserve curated mapping; repair only demonstrable negation handling.
            if 'non-neuron' in panel.lower():value='UNMAPPABLE'
            maps.append(dict(library=name,panel=panel,term=native_term(panel),L1=value,
                mapping_source='frozen mapping_L1_v3; non-neuron negation guard',
                mapping_source_sha256=sha(ROOT/'handoff/markers_v3/mapping_L1_v3.csv')))
    context=json.loads((OLD/'markers/brain_GBM.json').read_text())
    for name,panels in context.items():
        assert name not in libraries
        libraries[name]=panels
        sources.append(dict(library=name,source=str(OLD/'markers/brain_GBM.json'),
            sha256=sha(OLD/'markers/brain_GBM.json'),type='prespecified CellMarker2 tissue/context',
            context='Brain normal/disease; Blood, Blood vessel, Lymph node; related union; AllHuman',
            selection='Anatomy/context only; no prediction scores'))
        for panel in panels:
            term=native_term(panel)
            previous=prior_terms.get(term.strip().lower(),set())
            value=next(iter(previous)) if len(previous)==1 else semantic(term)
            if 'non-neuron' in term.lower():value='UNMAPPABLE'
            maps.append(dict(library=name,panel=panel,term=term,L1=value,
                mapping_source='same-term frozen mapping_L1_v3 when unique, otherwise prespecified semantic()',mapping_source_sha256=sha(__file__)))
    write_json(d/'libraries.json',libraries)
    pd.DataFrame(maps).to_csv(d/'panel_L1_mapping.csv',index=False)
    shutil.copy2(OLD/'markers/native_panel_metadata.csv',d/'CellMarker_native_panel_metadata.csv')
    meta=pd.read_csv(d/'CellMarker_native_panel_metadata.csv',keep_default_na=False)
    counts=[]
    for lib,panels in libraries.items():
        rows=[r for r in maps if r['library']==lib]
        supported=sorted({r['L1'] for r in rows}&set(L1))
        counts.append(dict(library=lib,n_panels=len(panels),n_unique_genes=len({g for v in panels.values() for g in v}),
            supported_L1=';'.join(supported),missing_L1=';'.join(sorted(set(L1)-set(supported))),
            n_unmapped_panels=sum(r['L1'] not in L1 for r in rows)))
    pd.DataFrame(counts).to_csv(d/'coverage_vocabulary.csv',index=False)
    src=ROOT/'handoff/refdb/Cell_marker_Human.xlsx'
    evidence=pd.read_excel(src).fillna('')
    # Keep original evidence metadata without pretending all markers are RNA-derived.
    evidence.to_csv(d/'CellMarker_source_records.csv.gz',index=False,compression='gzip')
    manifest=dict(status='completed',libraries=list(libraries),n_libraries=len(libraries),
        fixed_library='CM2_glioma_other',fixed_cutoff='mean',
        cutoff_grid=['none','mean','0.5'],native_labels_preserved_through_DL=True,
        mapping_applied_only_after_terminal_prediction=True,
        full_panel_gene_denominators_preserved=True,zero_overlap_panels_not_silently_removed=True,
        no_target_expression_or_label_based_marker_selection=True,
        evidence_columns=list(evidence.columns),sources=sources,
        statement='References may mix species and assay evidence. No heldout-patient independence is claimed for a database evidence source without source-study audit.',
        files={p.name:sha(p) for p in d.iterdir() if p.is_file()},
        source_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],completed_at=utc())
    write_json(d/'manifest.json',manifest);complete(d,'manifest.json','MARKERS_COMPLETE')
    print('MARKERS_COMPLETE',len(libraries),json.dumps(counts),flush=True)

if __name__=='__main__':run()
