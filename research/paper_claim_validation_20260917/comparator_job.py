"""Run frozen comparator implementations and score only after predictions finish."""
import json
import subprocess
import sys
from pathlib import Path
from common import OUT, RSCRIPT, require_slurm, checked, sha

def run(method,sample):
    require_slurm()
    import evaluate_comparator
    source=Path(__file__).resolve().parent
    if method=='scType':
        audit=OUT/'comparators/scType/TKU4163'
        assert checked(audit,'audit_sensitivity_manifest.json','AUDIT_SENSITIVITY_COMPLETE')
        am=json.loads((audit/'audit_sensitivity_manifest.json').read_text())
        assert am['source_sha256']==sha(source/'sctype_R.R')
        dest=OUT/'comparators/scType'/sample
        if not checked(dest,'cohort_manifest.json','COHORT_COMPLETE'):
            subprocess.run([RSCRIPT,str(source/'sctype_R.R'),sample,'cohort'],check=True)
    elif method=='SCINA':
        for i in range(16):
            assert checked(OUT/f'comparators/SCINA/TKU4163/L{i:02d}')
            if not checked(OUT/f'comparators/SCINA/{sample}/L{i:02d}'):
                subprocess.run([RSCRIPT,str(source/'scina_R.R'),sample,str(i)],check=True)
    elif method=='scCATCH':
        audit=OUT/'comparators/scCATCH/TKU4163/hvg2000/PCA30_SNN'
        assert checked(audit,'audit_manifest.json','AUDIT_COMPLETE')
        am=json.loads((audit/'audit_manifest.json').read_text())
        assert am['source_sha256']==sha(source/'sccatch_R.R')
        dest=OUT/'comparators/scCATCH'/sample/'hvg2000/UMAP2_HDBSCAN_R'
        if not checked(dest,'cohort_manifest.json','COHORT_COMPLETE'):
            subprocess.run([RSCRIPT,str(source/'sccatch_R.R'),sample,'cohort','hvg2000','UMAP2_HDBSCAN_R'],check=True)
    else:raise ValueError(method)
    evaluate_comparator.run(method,sample)

if __name__=='__main__':run(*sys.argv[1:])
