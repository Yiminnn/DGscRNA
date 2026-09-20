from pathlib import Path
import sys,os,json
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
sys.path.insert(0,str(ROOT/'handoff/paper_claim_validation_20260917'))
from common import OUT,L1,require_slurm,sha,write_json,checked,utc
DEST=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/comparison/scDeepSort_LogNormalize'
require_slurm()
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
parity=json.loads((DEST/'R_full_input_parity.json').read_text());assert parity['status']=='passed'
proof=[]
for sample in ['TKU4163','NL022']:
    dest=DEST/sample;assert checked(dest)
    manifest=json.loads((dest/'manifest.json').read_text())
    for name,digest in manifest['files'].items():assert sha(dest/name)==digest,name
    old=OUT/'comparators/scDeepSort'/sample
    oldman=json.loads((old/'manifest.json').read_text())
    assert sha(old/'predictions.csv.gz')==oldman['predictions_sha256']==manifest['old_predictions_untouched_sha256']
    inp=json.loads((dest/'input_manifest.json').read_text())
    assert inp['model_input_integers_cast'] is False and inp['normalized_exactly_once'] is True
    for path,digest in inp['pretrained_sources'].items():assert sha(path)==digest
    assert inp['pretrained_sources'][next(p for p in inp['pretrained_sources'] if p.endswith('human-Brain.pt'))]==oldman['model_sha256']
    assert parity['samples'][sample]['actual_predictor_csv_sha256']==inp['input_csv_sha256']
    labs=pd.read_csv(dest/'mapped_labels.csv.gz',dtype=str,keep_default_na=False)
    coords=pd.read_csv(OUT/'GBM'/sample/'hvg2000/UMAP2.csv',index_col=0)
    assert list(labs.cell_id)==list(coords.index)
    fig,axes=plt.subplots(1,3,figsize=(14,4.8),layout='constrained')
    colors=dict(zip(L1,['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9','#332288','#88CCEE','#44AA99','#AA4499','#999933']))
    colors.update(Unknown='#bbbbbb',UNMAPPABLE='#333333',AMBIGUOUS_NEURON='#777777',NO_L1_COUNTERPART='#555555')
    for ax,column,title in zip(axes,['L1','historical_raw_counts','corrected_LogNormalize'],['Author L1','Historical raw-count input','Corrected LogNormalize input']):
        for lab in sorted(set(labs[column])):
            mask=labs[column].eq(lab).to_numpy()
            ax.scatter(coords.iloc[mask,0],coords.iloc[mask,1],s=3,c=[colors.get(lab,'#777777')],linewidths=0,rasterized=True)
        ax.set_title(title);ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('Fixed display UMAP1');ax.set_ylabel('Fixed display UMAP2')
    present=set(labs[['L1','historical_raw_counts','corrected_LogNormalize']].to_numpy().ravel())
    fig.legend(handles=[Line2D([],[],marker='o',ls='',color=colors.get(lab,'#777777'),label=lab,markersize=5) for lab in L1+sorted(present-set(L1)) if lab in present],
        loc='outside lower center',ncol=6,frameon=False,fontsize=8)
    fig.suptitle(sample+': same pretrained scDeepSort weights and native threshold\nInput-normalization repair; published Brain vocabulary still lacks malignant labels',fontsize=12)
    for ext in ['png','pdf','svg']:fig.savefig(DEST/f'{sample}_input_comparison.{ext}',dpi=260,bbox_inches='tight')
    plt.close(fig)
    proof.append(dict(sample=sample,n_cells=len(labs),old_predictions_unchanged=True,
        identical_published_checkpoint=True,every_CSV_value_matches_R=True))
write_json(DEST/'verification.json',dict(status='passed',job=os.environ['SLURM_JOB_ID'],completed_at=utc(),samples=proof,
    source_sha256=sha(__file__),R_parity_sha256=sha(DEST/'R_full_input_parity.json')))
state=json.loads((DEST/'status.json').read_text());state.update(verified=True,updated_at=utc())
state['jobs'].append(os.environ['SLURM_JOB_ID']);state['evidence'].extend([str(DEST/'verification.json'),str(DEST/'R_full_input_parity.json')])
write_json(DEST/'status.json',state)
print(json.dumps(proof,indent=2))
