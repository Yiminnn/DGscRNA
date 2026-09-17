"""Check the small-patient inference contract and presentation dependencies."""
import os
from pathlib import Path
from ptc_common import BASE,require_slurm,sha,utc,write_json
require_slurm()
import numpy as np
import pandas as pd
import matplotlib
import nbformat
import nbclient
import tabulate
from summarize_ptc import effect_values,paired,holm

a=effect_values([1,1,1,1]);assert a['p_exact']==.125 and a['delta_mean']==a['ci_lower']==a['ci_upper']==1
a=effect_values([1,-1,1,-1]);assert a['p_exact']==1 and a['delta_mean']==0
q=pd.DataFrame({'sample':list('abcdefgh'),'patient':['P1','P1','P2','P2','P3','P3','P4','P4'],
                'score':[0,2,0,4,0,6,0,8]})
b=q.copy();b['score']=0
e=paired(q,b,'check','check',metrics=['score'])[0]
assert e['n_paired_samples']==8 and e['n_paired_patients']==4 and e['delta_mean']==2.5
expected=np.quantile(np.array([np.mean(x) for x in __import__('itertools').product([1.,2.,3.,4.],repeat=4)]),[.025,.975])
np.testing.assert_array_equal([e['ci_lower'],e['ci_upper']],expected)
h=holm(pd.DataFrame({'family':['x']*5,'p_exact':[.125,.25,.5,.75,1.]}),['family'])
np.testing.assert_allclose(h.p_holm,[.625,1,1,1,1])
out=BASE/'verification/analysis_preflight';out.mkdir(exist_ok=True)
write_json(out/'manifest.json',dict(status='passed',checks=['four-patient minimum p=.125','symmetric zero effect',
    'sample differences averaged within patient','exhaustive bootstrap quantiles','five-comparison Holm','presentation dependencies'],
    source_sha256=sha(Path(__file__)),statistics_source_sha256=sha(Path(__file__).with_name('summarize_ptc.py')),
    job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
(out/'COMPLETE').write_text(sha(out/'manifest.json')+'\n')
print('Patient-level inference and notebook/figure dependencies checked',flush=True)
