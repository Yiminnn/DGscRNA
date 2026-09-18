"""Preserve shared marker genes; recover gene-specific evidence IDs where available."""
import json
import os
import re
from common import OUT, require_slurm, write_json, sha, complete, utc

def run():
    require_slurm()
    import pandas as pd
    root=OUT/'markers/scCATCH';root.mkdir(exist_ok=True)
    libs=json.loads((OUT/'markers/libraries.json').read_text())
    t=pd.read_csv(OUT/'markers/CellMarker_source_records.csv.gz',dtype=str,keep_default_na=False)
    for c in t:t[c]=t[c].str.strip()
    t=t[t.species.str.lower().eq('human')].copy()
    t['native']='CM2+'+t.cell_type+'+'+t.tissue_class+'+'+t.cancer_type.replace('','unspecified')+'+'+t.cell_name
    exact={};collapsed={}
    for r in t.itertuples():
        source=str(r.PMID).removesuffix('.0') or 'CellMarker-record-without-PMID'
        contexts=['all']
        if r.cell_type=='Normal cell':contexts+=['CM2_brain_normal']
        if 'glioblast' in r.cancer_type.lower():contexts+=['CM2_glioblastoma']
        if 'glioma' in r.cancer_type.lower() and 'glioblast' not in r.cancer_type.lower():contexts+=['CM2_glioma_other']
        for gene in re.split('[,;]',r.Symbol):
            gene=gene.strip()
            if not gene:continue
            exact.setdefault((r.native,gene),set()).add(source)
            for context in contexts:collapsed.setdefault((context,r.tissue_class,r.cell_name,gene),set()).add(source)
    manifests=[]
    for i,(lib,panels) in enumerate(libs.items()):
        rows=[];matched=0;fallback=0
        for panel,genes in panels.items():
            for gene in genes:
                pmids=exact.get((panel,gene),set())
                if not pmids and panel.startswith('CM2_'):
                    parts=panel.split('+',2)
                    if len(parts)==3:pmids=collapsed.get((parts[0],parts[1],parts[2],gene),set())
                if pmids:matched+=1
                else:
                    # This is explicitly a source-library evidence unit, not a PMID
                    # or an invented count of supporting research articles.
                    pmids={'library-source:'+panel.split('+',1)[0]};fallback+=1
                for evidence in sorted(pmids):
                    rows.append(dict(gene=gene,celltype=panel,pmid=evidence,
                        subtype1=None,subtype2=None,subtype3=None))
        p=root/f'L{i:02d}.csv.gz';pd.DataFrame(rows).drop_duplicates().to_csv(p,index=False,compression='gzip')
        manifests.append(dict(library=lib,path=str(p),sha256=sha(p),gene_panel_pairs_with_gene_specific_evidence=matched,
            gene_panel_pairs_with_library_source_unit=fallback,genes_and_panels_identical_to_frozen_DG_library=True))
    write_json(root/'manifest.json',dict(status='completed',libraries=manifests,
        evidence='Native CellMarker panels use exact panel-gene PMID rows. Curated CellMarker panels match their named tissue/context and genes. Other unrecovered records use a single explicitly named library-source unit, not fabricated PMIDs.',
        limitation='Library-level sources do not recover original gene-level article evidence; the matched-panel benchmark must retain this limitation.',
        input_sha256=sha(OUT/'markers/libraries.json'),source_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(root)

if __name__=='__main__':run()
