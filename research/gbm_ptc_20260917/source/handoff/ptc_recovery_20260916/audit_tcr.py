#!/usr/bin/env python3
"""Audit TCR detection definitions against original contigs and supplied S2.

Detection-negative means no evidence under that rule, not a known non-T cell.
Rules are source/quality definitions fixed before comparing annotation labels.
"""
from pathlib import Path
import sys
import time
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'hvg_ptc_20260916'))
from common import ROOT, OUT, require_slurm, sha, utc, write_json, runtime_record


def run():
    require_slurm()
    import pandas as pd
    base = OUT / 'ptc_recovery'
    assert (base / 'ARCHIVE_STAGED').exists()
    archive = base / 'archive/tcr'
    dest = base / 'tcr_audit'
    dest.mkdir(exist_ok=True)
    start = time.monotonic()
    meta = pd.read_csv(archive / 'rawdata/metadata.txt', sep='\t').set_index('Sample')
    s2path = ROOT / 'paper/submission_v16/Supplementary table S2_T cell type identification_DGscRNA.xlsx'
    s2 = pd.read_excel(s2path, header=3)
    s2.columns = s2.columns.astype(str).str.strip()
    s2 = s2.loc[s2.Sample_ID.notna()].set_index('Sample_ID')
    assert s2.index.is_unique
    s2truth = pd.to_numeric(s2.is_T_cell_Real, errors='raise').astype(int)
    provenance, summaries, allcells = [], [], []
    def true(series):
        return series.astype(str).str.lower().isin(['true','t','1'])
    for sample in meta.index:
        path = archive / f'rawdata/TCR/{sample}/filtered_contig_annotations.csv'
        d = pd.read_csv(path)
        for col in ['barcode','chain','productive','high_confidence','is_cell']:
            assert col in d, (sample,col,d.columns.tolist())
        d['canonical_cell_id'] = sample + '_' + d.barcode.astype(str).str.split('-').str[0]
        knownchain = d.chain.isin(['TRA','TRB','TRG','TRD'])
        productive = true(d.productive) & knownchain
        high = productive & true(d.high_confidence)
        eligible = high & true(d.is_cell)
        masks = dict(any_filtered_contig=pd.Series(True,index=d.index),
            any_TCR_chain=knownchain, any_productive_TCR=productive,
            high_confidence_productive_TCR=high,
            cell_high_confidence_productive_TCR=eligible,
            productive_TRB=productive & d.chain.eq('TRB'))
        cells = pd.DataFrame(index=pd.Index(d.canonical_cell_id.unique(), name='cell_id'))
        cells['sample'] = sample
        cells['patient'] = meta.loc[sample,'Patient']
        cells['n_contigs'] = d.groupby('canonical_cell_id').size().reindex(cells.index)
        for name, mask in masks.items():
            cells[name] = cells.index.isin(d.loc[mask,'canonical_cell_id'])
        a = set(d.loc[eligible & d.chain.eq('TRA'),'canonical_cell_id'])
        b = set(d.loc[eligible & d.chain.eq('TRB'),'canonical_cell_id'])
        cells['paired_productive_TRA_TRB'] = cells.index.isin(a & b)
        rules = list(masks) + ['paired_productive_TRA_TRB']
        target_ids = s2.index[s2.index.str.startswith(sample+'_')]
        target = pd.DataFrame(index=target_ids)
        target['sample'] = sample
        target['patient'] = meta.loc[sample,'Patient']
        target['S2_is_T_cell_Real'] = s2truth.loc[target_ids]
        for rule in rules:
            target[rule] = cells[rule].reindex(target_ids, fill_value=False).astype(bool)
            agree = target[rule].astype(int).eq(target.S2_is_T_cell_Real)
            summaries.append(dict(sample=sample, patient=meta.loc[sample,'Patient'],
                rule=rule, n_contig_rows=len(d), n_TCR_barcodes=len(cells),
                n_rule_positive_all_TCR=int(cells[rule].sum()), n_S2_cells=len(target),
                n_S2_positive=int(target.S2_is_T_cell_Real.sum()),
                n_rule_positive_S2=int(target[rule].sum()),
                n_exact_S2=int(agree.sum()), n_disagree_S2=int((~agree).sum()),
                n_S2_positive_rule_negative=int((target.S2_is_T_cell_Real.eq(1)&~target[rule]).sum()),
                n_S2_negative_rule_positive=int((target.S2_is_T_cell_Real.eq(0)&target[rule]).sum())))
        target.to_csv(dest / f'{sample}.S2_detection_comparison.csv.gz')
        cells.to_csv(dest / f'{sample}.original_TCR_barcode_rules.csv.gz')
        allcells.append(target)
        provenance.append(dict(sample=sample, source=str(path), sha256=sha(path),
            columns=list(d.columns), n_contigs=len(d), n_unique_barcodes=len(cells),
            is_cell_values=d.is_cell.astype(str).value_counts().to_dict(),
            productive_values=d.productive.astype(str).value_counts().to_dict(),
            high_confidence_values=d.high_confidence.astype(str).value_counts().to_dict()))
        print(sample,'TCR barcode and S2 comparison finished',flush=True)
    pd.concat(allcells).loc[s2.index].to_csv(dest / 'all_S2_TCR_detection_rules.csv.gz')
    summary = pd.DataFrame(summaries)
    summary.to_csv(dest / 'detection_rule_comparison_by_sample.csv', index=False)
    countcols=[c for c in summary if c.startswith('n_')]
    summary.groupby('rule')[countcols].sum().to_csv(dest / 'detection_rule_comparison_cohort.csv')
    write_json(dest / 'manifest.json', dict(status='completed',timestamp=utc(),
        n_S2=len(s2), S2_sha256=sha(s2path), original_files=provenance,
        primary_quality_rule='cell_high_confidence_productive_TCR',
        rule_interpretation='Detection evidence; undetected is not a known non-T cell.',
        S3_status='separate reconciliation required',
        source_sha256=sha(Path(__file__)), seconds=time.monotonic()-start, **runtime_record()))
    (dest / 'COMPLETE').write_text(sha(dest / 'manifest.json')+'\n')


if __name__ == '__main__':
    run()
