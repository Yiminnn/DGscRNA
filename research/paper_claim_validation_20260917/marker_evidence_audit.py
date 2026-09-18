"""Expose original source/assay records and unresolved provenance without relabeling."""
import json
import os
from common import ROOT, OUT, require_slurm, checked, complete, sha, write_json, utc

def run():
    require_slurm()
    import pandas as pd
    target=OUT/'marker_evidence_summary';target.mkdir(exist_ok=True)
    libraries=json.loads((OUT/'markers/libraries.json').read_text())
    evidence=pd.read_csv(OUT/'markers/CellMarker_source_records.csv.gz',dtype=str,keep_default_na=False)
    evidence['PMID']=evidence.PMID.str.strip().str.removesuffix('.0')
    gene_records=evidence.assign(gene=evidence.Symbol.str.split('[,;]',regex=True)).explode('gene')
    gene_records['gene']=gene_records.gene.str.strip()
    metadata=gene_records[['PMID','gene','species','tissue_class','tissue_type','cancer_type','technology_seq','marker_source','Title','journal','year']].drop_duplicates()
    coverage=pd.read_csv(OUT/'markers/coverage_vocabulary.csv')
    detail=[];summary=[]
    for i,lib in enumerate(libraries):
        source=pd.read_csv(OUT/'markers/scCATCH'/f'L{i:02d}.csv.gz',dtype=str,keep_default_na=False)
        pairs=source[['celltype','gene','pmid']].drop_duplicates().rename(columns={'pmid':'PMID'})
        paired=pairs.merge(metadata,on=['PMID','gene'],how='left')
        paired.insert(0,'library',lib)
        paired['gene_specific_source_recovered']=paired.Title.notna()
        detail.append(paired)
        found=paired[paired.gene_specific_source_recovered]
        missing=pairs[pairs.PMID.str.startswith('library-source:')]
        # PMID absence rules out this exact indexed article only. It cannot rule
        # out participant overlap, an earlier article, or shared author signatures.
        target_hits=sorted(set(found.PMID)&{'40346361','40346362'})
        original_overlap=lib in ['CARE_TME','BrainAtlas112','UNION_all']
        summary.append(dict(library=lib,n_panels=len(libraries[lib]),
            n_gene_panel_pairs=len(pairs[['celltype','gene']].drop_duplicates()),
            n_pairs_with_gene_specific_evidence=len(found[['celltype','gene']].drop_duplicates()),
            n_pairs_with_library_only_evidence=len(missing[['celltype','gene']].drop_duplicates()),
            n_indexed_PMIDs=found.PMID.nunique(),species=';'.join(sorted(set(found.species.dropna())-{''})),
            technology_seq=';'.join(sorted(set(found.technology_seq.dropna())-{''})),
            evidence_methods=';'.join(sorted(set(found.marker_source.dropna())-{''})),
            target_article_PMIDs=';'.join(target_hits),author_label_construction_overlap=original_overlap,
            interpretation='Author-label concordance; excluded from primary marker selection' if original_overlap else
                'No exact target-article hit' if not target_hits else 'Target-article overlap requires concordance interpretation',
            patient_overlap_verified=False,all_markers_RNA_derived_verified=False))
    table=pd.DataFrame(summary).merge(coverage.drop(columns='n_panels'),on='library',validate='one_to_one')
    table.to_csv(target/'library_source_assay_coverage_audit.csv',index=False)
    pd.concat(detail,ignore_index=True).to_csv(target/'gene_panel_source_records.csv.gz',index=False,compression='gzip')
    columns=['library','PMID','Title','journal','year','species','technology_seq','marker_source']
    pd.concat(detail,ignore_index=True)[columns].drop_duplicates().to_csv(target/'source_studies.csv',index=False)
    report='''# Marker provenance and interpretation

The frozen libraries preserve their original panel names and genes. The attached tables recover gene-specific CellMarker article records where the archived database supports that link; missing links remain library-level evidence. Mixed or blank assay metadata are not recoded as RNA evidence. Species and assay summaries describe recovered records only, not an independence guarantee for every gene.

CARE_TME and BrainAtlas112 participated in the original author annotation process; UNION_all contains them. All three remain visible as concordance results and are excluded from primary marker selection. The archived source report is `handoff/FROZEN_single_sample_TKU3186/report/REPORT_TKU3186_for_advisor.md`, sections 4–5. Other signatures from that repository are external-reference candidates, but a full per-participant independence audit has not been established.

The two GSE274546 publication PMIDs are 40346361 and 40346362. Their absence from an archived database table only excludes those exact indexed articles; it does not establish absence of shared participants or earlier reused data. The primary GBM paper identifies the dataset and author analysis code: [GSE274546 publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC12081307/), [author repository](https://github.com/dravishays/GBM-CARE-WT).

This audit does not replace the frozen marker roster after viewing its accuracy. Missing cell classes remain errors in the primary all-cell metric. A normal/developmental brain reference is not assumed to cover malignant or adult immune cells merely because its tissue name contains brain.
'''
    (target/'MARKER_EVIDENCE.md').write_text(report)
    write_json(target/'manifest.json',dict(status='completed',n_libraries=16,
        source_file_sha256=sha(OUT/'markers/CellMarker_source_records.csv.gz'),
        original_source_report_sha256=sha(ROOT/'handoff/FROZEN_single_sample_TKU3186/report/REPORT_TKU3186_for_advisor.md'),
        source_urls=['https://pmc.ncbi.nlm.nih.gov/articles/PMC12081307/','https://github.com/dravishays/GBM-CARE-WT'],
        independence_limitations_retained=True,files={p.name:sha(p) for p in target.iterdir() if p.suffix in ['.csv','.gz','.md']},
        job=os.environ['SLURM_JOB_ID'],source_sha256=sha(__file__),completed_at=utc()))
    complete(target)

if __name__=='__main__':run()
