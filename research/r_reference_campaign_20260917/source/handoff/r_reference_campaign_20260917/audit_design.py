"""Compare actual gene universes and cohort cell-order contracts, without relabeling."""
import os,json,hashlib
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import pandas as pd
    rows=[]
    for group in ['NMT','TTU']:
        parent=OUT/'PTC_ablation'
        comparisons=[(f'PTC_{group}_CCA2000',f'PTC_{group}_NONE_RNA2000','CCA versus RNA'),
                     (f'PTC_{group}_NONE_RNA2000',f'PTC_{group}_HARMONY_RNA2000','RNA versus Harmony')]
        comparisons += [(f'PTC_{group}_GEOMETRY2000_FIXED_CCAall_DL2000',f'PTC_{group}_GEOMETRY{h}_FIXED_CCAall_DL2000','isolated geometry') for h in ['500','1000','3000','5000','all']]
        for left,right,meaning in comparisons:
            a,b=parent/left,parent/right
            if not (a/'PREPARED').exists() or not (b/'PREPARED').exists():continue
            ga=(a/'DL_features.txt').read_text().splitlines();gb=(b/'DL_features.txt').read_text().splitlines()
            ma=json.loads((a/'prepare_manifest.json').read_text());mb=json.loads((b/'prepare_manifest.json').read_text())
            ca=pd.read_csv(a/'cells.csv',dtype=str).cell_id;cb=pd.read_csv(b/'cells.csv',dtype=str).cell_id
            same_binary=ma['DL_binary_sha256']==mb['DL_binary_sha256']
            same_order=ca.tolist()==cb.tolist()
            r=dict(group=group,comparison=meaning,left=left,right=right,
              left_genes=len(ga),right_genes=len(gb),overlap=len(set(ga)&set(gb)),
              same_gene_set=set(ga)==set(gb),same_gene_order=ga==gb,same_cell_order=same_order,
              identical_DL_binary_sha256=same_binary)
            if meaning in ['RNA versus Harmony','isolated geometry']:
                assert same_binary and same_order and ga==gb,r
            rows.append(r)
    dest=OUT/'verification';dest.mkdir(exist_ok=True)
    pd.DataFrame(rows).to_csv(dest/'actual_feature_and_input_contrasts.csv',index=False)
    report=dict(job=os.environ['SLURM_JOB_ID'],status='passed_for_available_preparations',comparisons=rows,
      caution='CCA versus RNA is a combined expression/geometry intervention; if actual feature sets differ this also changes gene identity. RNA/Harmony and isolated-geometry contrasts require byte-identical DL input.')
    (dest/'design_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':run()
