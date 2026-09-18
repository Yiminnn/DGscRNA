"""Isolated, byte-identical input aliases for independent resource repetitions."""
import json
import os
from common import OUT,require_slurm,sha,write_json,complete,checked

def output(method,n,repeat=0):
    p=OUT/'scalability'/str(n) if repeat==0 else OUT/'scalability_repeats'/f'repeat{repeat}'/str(n)
    return p if method=='DG-scRNA' else p/method

def input_alias(method,n,repeat=0):
    require_slurm()
    if repeat==0:return f'SCALE_{n}'
    assert repeat in [1,2]
    original=OUT/'inputs'/f'SCALE_{n}'
    assert checked(original,'input_manifest.json','INPUT_COMPLETE')
    name=f'SCALE_{n}_{method.replace("-","")}_repeat{repeat}'
    dest=OUT/'inputs'/name;dest.mkdir(parents=True,exist_ok=True)
    for file in ['x.bin','i.bin','p.bin','genes.csv','cells_fit.csv']:
        target=dest/file
        if not target.exists():os.link(original/file,target)
        assert os.stat(target).st_ino==os.stat(original/file).st_ino
    m=json.loads((original/'input_manifest.json').read_text())
    m.update(sample=name,resource_repeat=repeat,resource_method=method,
        input_alias_of=str(original),input_alias_manifest_sha256=sha(original/'input_manifest.json'),
        count_gene_cell_files_hardlinked_byte_identical=True)
    if checked(dest,'input_manifest.json','INPUT_COMPLETE'):
        assert json.loads((dest/'input_manifest.json').read_text())==m
    else:write_json(dest/'input_manifest.json',m);complete(dest,'input_manifest.json','INPUT_COMPLETE')
    return name
