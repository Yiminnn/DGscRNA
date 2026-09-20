"""Independent accounting checks on the saved D/F supplement (SLURM only)."""
from pathlib import Path
import os
import json
import sys

ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
sys.path.insert(0,str(ROOT/'handoff/paper_claim_validation_20260917'))
from common import require_slurm,sha,checked,write_json,utc,L1
OUT=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/evidence'

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    assert checked(OUT)
    manifest=json.loads((OUT/'manifest.json').read_text())
    for name,digest in manifest['files'].items():assert sha(OUT/name)==digest,name
    for name,digest in manifest['sources'].items():assert sha(ROOT/name)==digest,name
    pairs=pd.read_csv(OUT/'pairwise_agreement_sample.csv.gz')
    expected={'all121':(121,59),'primary97':(97,55)}
    checks=[]
    for cohort,(samples,patients) in expected.items():
        p=pairs[pairs.cohort==cohort]
        assert p['sample'].nunique()==samples and p.patient.nunique()==patients
        assert len(p)==samples*15*3
        assert not p.duplicated(['sample','method_a','method_b','subset']).any()
        assert (p.n_subset<=p.n_all_cells).all()
        assert (p.n_same_supported_label<=p.n_same_label).all()
        assert (p.n_same_label<=p.n_subset).all()
        cells=0
        for f in (OUT/'predictions'/cohort).glob('*.csv.gz'):
            frame=pd.read_csv(f)
            assert frame.cell_id.is_unique
            for group,methods in [('all_six',['DG-scRNA','scType','scCATCH','SCINA','SingleR','scDeepSort']),
                                  ('four_marker_methods',['DG-scRNA','scType','scCATCH','SCINA'])]:
                a=frame[methods].to_numpy()
                expected_unanimous=(a==a[:,:1]).all(1)&np.isin(a,L1).all(1)
                np.testing.assert_array_equal(frame[group+'_unanimous_supported'],expected_unanimous)
                count=np.stack([(a==label).sum(1) for label in L1],axis=1)
                assert ((frame[group+'_majority']!='No_majority').to_numpy()==(count.max(1)>len(methods)/2)).all()
            cells+=len(frame)
        assert len(list((OUT/'predictions'/cohort).glob('*.csv.gz')))==samples
        checks.append(dict(cohort=cohort,samples=samples,patients=patients,n_cells=cells,pair_rows=len(p)))
    g=pd.read_csv(OUT/'consensus_strata_sample.csv')
    wide=g.pivot(index=['cohort','sample','group'],columns='category',values='n_subset')
    np.testing.assert_array_equal(wide.all_cells,wide.unanimous_supported+wide.nonunanimous_or_unsupported)
    np.testing.assert_array_equal(wide.nonunanimous_or_unsupported,wide.any_unsupported+wide.all_supported_but_disagree)
    provenance=json.loads((OUT/'NL022_expression_provenance.json').read_text())
    assert provenance['normalization_checked_cells']==3628 and provenance['normalization_checked_genes']==2000
    assert provenance['normalization_all_HVG_max_abs_R_difference']<1e-6
    values=pd.read_csv(OUT/'NL022_display_gene_cell_expression.csv.gz')
    labels=pd.read_csv(OUT/'NL022_frozen_endpoint_labels.csv.gz')
    assert list(values.cell_id)==list(labels.cell_id)
    assert all((OUT/f'{stem}.{ext}').stat().st_size>1000
        for stem in ['GBM_cross_tool_agreement','GBM_consensus_vs_author_L1','NL022_selected_marker_dotplot',
                     'NL022_selected_marker_violin','NL022_markers_initial_terminal_same_coordinates']
        for ext in ['png','pdf','svg'])
    note=dict(status='passed',verified_at=utc(),job=os.environ['SLURM_JOB_ID'],cohorts=checks,
        files_and_sources_sha256_verified=True,all_consensus_cell_rows_rechecked=True,
        independent_exhaustive_denominator_partitions=True,NL022_R_normalization_parity=True,
        source_sha256=sha(__file__))
    write_json(OUT/'verification.json',note)
    print(json.dumps(note,indent=2))

if __name__=='__main__':run()
