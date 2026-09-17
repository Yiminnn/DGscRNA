"""Package derived results in the already authorized OneDrive directory layout."""
import ast,fcntl,gzip,hashlib,json,os,shutil,subprocess,tarfile
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'
CODE=ROOT/'handoff/r_reference_campaign_20260917'
PARENT=ROOT/'results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916'
STAGE=PARENT/'GBM_PTC_results_20260916'
REMOTE='onedrive:work_od/share/dgscrna_GSE274546_TKU3186/01_report/GBM_PTC_results_20260916'

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import pandas as pd
    summary=json.loads((OUT/'summary/campaign_summary.json').read_text());assert summary['status']=='complete'
    nb=json.loads((OUT/'summary/notebook_update_manifest.json').read_text());assert nb['status']=='complete'
    assert sha(ROOT/'notebooks/dgscrna_results.ipynb')==nb['current_sha256']
    inv=pd.read_csv(OUT/'summary/completion_inventory.csv');assert len(inv)==summary['expected_analysis_units'] and inv.audited.all()
    design=json.loads((OUT/'verification/design_audit.json').read_text());assert len(design['comparisons'])==14
    early=list((OUT/'verification/pre_guard_density').glob('*.json'));assert len(early)==16
    assert all(json.loads(p.read_text())['status']=='passed' for p in early)
    seeds=json.loads((OUT/'verification/PTC_marker_seed_coverage_audit.json').read_text())
    assert seeds['status']=='complete' and seeds['selected_group_routes']==32
    assert all(json.loads((OUT/'verification/ptc_density'/f'PTC_{g}_CCA2000.json').read_text())['status']=='passed' for g in ['NMT','TTU'])
    batch_dir=OUT/'verification/batch_biology'
    batch=json.loads((batch_dir/'manifest.json').read_text())
    assert batch['status']=='complete' and batch['n_units']==6 and batch['n_unit_spaces']==12
    assert batch['same_queries_and_cell_order_across_correction_arms']
    assert (batch_dir/'COMPLETE').read_text().strip()==sha(batch_dir/'manifest.json')
    assert all(sha(batch_dir/name)==value for name,value in batch['outputs'].items())
    task_recovery=json.loads((OUT/'verification/terminal_tasklist_recovery.json').read_text())
    assert task_recovery['status']=='passed' and task_recovery['terminal_arms_reaudited']==204
    grid=json.loads((OUT/'verification/frozen_grid_audit.json').read_text())
    assert grid['status']=='passed' and grid['units']==100
    score_recovery=json.loads((OUT/'verification/score_recovery_and_DEG_workers.json').read_text())
    assert score_recovery['status']=='passed'
    accounting=json.loads((OUT/'resources/accounting_manifest.json').read_text())
    assert accounting['status']=='collected' and accounting['job']==os.environ['SLURM_JOB_ID']
    for name,value in accounting['files'].items():assert sha(OUT/'resources'/name)==value
    for entry in accounting['failed_logs']:assert sha(OUT/entry['path'])==entry['sha256']
    code_files=sorted(p for p in CODE.rglob('*') if p.is_file() and '__pycache__' not in str(p) and p.name!='CURRENT_TASK.md' and not p.name.endswith(('.pyc','.next')))
    for p in code_files:
        if p.suffix=='.py':ast.parse(p.read_text(),filename=str(p))
        if p.suffix=='.sbatch':subprocess.run(['bash','-n',str(p)],check=True)
    rs=[str(p) for p in code_files if p.suffix=='.R']
    subprocess.run(['/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript','-e',
       'for(p in commandArgs(trailingOnly=TRUE))parse(p)',*rs],check=True)
    checked=dict(status='passed',python_files=sum(p.suffix=='.py' for p in code_files),
      R_files=len(rs),shell_files=sum(p.suffix=='.sbatch' for p in code_files),job=os.environ['SLURM_JOB_ID'])
    (OUT/'verification/source_syntax_checks.json').write_text(json.dumps(checked,indent=2)+'\n')
    source_index=[dict(path=str(p.relative_to(ROOT)),sha256=sha(p)) for p in code_files]
    (OUT/'summary/source_manifest.json').write_text(json.dumps(source_index,indent=2)+'\n')
    files=[];artifacts=[];cached_deg_aliases=[]
    def copy(p):
        rel=p.relative_to(ROOT);target=STAGE/rel;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(p,target);files.append(str(rel))
    def bundle(paths,target,relative_root,extra_members=()):
        target.parent.mkdir(parents=True,exist_ok=True)
        temp=target.with_name(target.name+'.part')
        with tarfile.open(temp,'w:gz',compresslevel=1) as archive:
            for p in sorted(paths):archive.add(p,arcname=str(p.relative_to(relative_root)),recursive=False)
            for source,archive_name in extra_members:archive.add(source,arcname=archive_name,recursive=False)
        temp.replace(target);files.append(str(target.relative_to(STAGE)))
    copy(ROOT/'notebooks/dgscrna_results.ipynb')
    for p in code_files:copy(p)
    for rel in ['handoff/ptc_recovery_20260916/refine.py','handoff/ptc_recovery_20260916/ptc_common.py',
                'handoff/ptc_recovery_20260916/label_rules.py','handoff/harmonize.py','handoff/deck_datasets_provenance.md']:
        copy(ROOT/rel)
    for sub in ['summary','markers','verification']:
        for p in sorted((OUT/sub).rglob('*')):
            if (p.is_file() and p.name not in ['DELIVERY_RECEIPT.json','DELIVERY_MANIFEST.json']
                and not p.name.startswith('REMOTE_RECEIPT_UPLOADED')):copy(p)
    # Full task/step resource accounting and failed-attempt logs are delivery artifacts.
    bundle([p for p in (OUT/'resources').rglob('*') if p.is_file() and not p.name.endswith('.part')],
      STAGE/OUT.relative_to(ROOT)/'slurm_accounting_and_failed_logs.tar.gz',OUT/'resources')
    for name in ['accounting_manifest.json','slurm_job_ledger.csv','failed_job_logs.json']:
        copy(OUT/'resources'/name)
    for p in OUT.glob('*.json'):
        if p.name not in ['dispatch_state.json']:copy(p)
    for p in OUT.glob('*.txt'):copy(p)
    for p in OUT.glob('*events.jsonl'):copy(p)
    for dataset_dir in sorted((OUT/'inputs').iterdir()):
        if dataset_dir.is_dir():
            selected=[p for p in dataset_dir.rglob('*') if p.is_file() and p.suffix not in ['.bin']]
            target=STAGE/dataset_dir.relative_to(ROOT)/'evaluation_inputs_and_manifests.tar.gz'
            bundle(selected,target,dataset_dir)
    for row in inv.itertuples(index=False):
        prep=Path(row.directory)
        for checkpoint in prep.glob('*/DEG_checkpoints/ASSEMBLED.json'):
            assembled=json.loads(checkpoint.read_text())
            validation_path=OUT/'verification/checkpointed_DEG_equivalence.json'
            validation=json.loads(validation_path.read_text())
            assert validation['status']=='passed' and validation['exact_all_fields_attributes_and_row_order']
            assert assembled['status']=='assembled' and assembled['validation_sha256']==sha(validation_path)
            assert assembled['source_sha256']==validation['source_sha256']
            assert assembled['installed_FindAllMarkers_sha256']==validation['installed_FindAllMarkers_sha256']
            assert assembled['sha256']==sha(checkpoint.parent.parent/'DEG.rds')
        paths=[]
        for p in prep.rglob('*'):
            if not p.is_file():continue
            relative=p.relative_to(prep);parts=relative.parts
            # Dense expression/anchor inputs and model weights remain on HPC.
            if p.name in ['expression_PCA30.rds','RNA_normalized.rds','anchors.rds','DL.float32.bin','model_state.pt','terminal.npz']:
                continue
            if (p.name.endswith(('.native_order_backup','.part')) or '.part.' in p.name
                or '.invalid_' in p.name):continue
            if 'figures' in parts or 'verification' in parts:
                copy(p);continue
            paths.append(p)
            if len(parts)==1 or (parts[0]=='evaluation' and len(parts)==2):copy(p)
        target=STAGE/prep.relative_to(ROOT)/'terminal_results_and_audits.tar.gz'
        # The two validated archived PTC routes reuse DEG files outside their unit
        # directory. Include those exact cached files in the delivery as well.
        extra_members=[]
        for manifest in sorted(prep.glob('*/score_manifest.json')):
            score=json.loads(manifest.read_text())
            if not score.get('DEG_file'):continue
            source=Path(score['DEG_file'])
            if source.is_relative_to(prep):continue
            assert source.is_file() and sha(source)==score['DEG_sha256']
            archive_name=str(manifest.parent.relative_to(prep)/'DEG_full_integrated2000.rds')
            existing=prep/archive_name
            if existing.exists():assert sha(existing)==score['DEG_sha256']
            else:extra_members.append((source,archive_name))
            cached_deg_aliases.append(dict(unit=row.unit,source=str(source),archived_as=archive_name,sha256=score['DEG_sha256']))
        bundle(paths,target,prep,extra_members)
        for p in prep.glob('*/terminal/*/training_manifest.json'):
            m=json.loads(p.read_text())
            for name in ['terminal.npz','model_state.pt']:
                if name in m['outputs']:artifacts.append(dict(path=str(p.parent/name),sha256=m['outputs'][name],
                    availability='retained on HPC; terminal label CSV and manifest included in per-unit delivery archive'))
        print('PACKAGED',row.unit,flush=True)
    # Executed per-job sources are compactly archived instead of thousands of loose files.
    bundle([p for p in (OUT/'execution_sources').rglob('*') if p.is_file()],
      STAGE/OUT.relative_to(ROOT)/'executed_sources.tar.gz',OUT/'execution_sources')
    local_record=OUT/'summary/local_retained_model_artifacts.jsonl.gz'
    with gzip.open(local_record,'wt') as f:
        for r in artifacts:f.write(json.dumps(r)+'\n')
    copy(local_record)
    deg_record=OUT/'summary/cached_DEG_delivery_sources.json'
    deg_record.write_text(json.dumps(cached_deg_aliases,indent=2)+'\n');copy(deg_record)
    layout=OUT/'summary/DELIVERY_LAYOUT.md'
    layout.write_text('''# Delivery layout

The existing `notebooks/dgscrna_results.ipynb` includes all prior cells plus the executed campaign section.
`summary/` contains the full grid, paired ablation differences, marker roster, label mappings and interpretation.
Each analysis unit has directly viewable figures, aggregate evaluation tables and an artifact-audit report.
Its `terminal_results_and_audits.tar.gz` contains all initial/final per-cell annotation CSVs, confidence values,
training histories/manifests, cluster assignments, DEG/density outputs and per-condition confusion matrices.
Extract that archive inside the unit directory to restore the detailed result tree.

Models and probability NPZ files remain on HPC, with original paths/checksums in
`local_retained_model_artifacts.jsonl.gz`; exported final labels for every arm are included in the archives.
Dense expression/anchor matrices and raw h5ad files are not duplicated into OneDrive.
Input sources and checksums remain in per-cohort input-manifest archives. Executed source copies are in
`executed_sources.tar.gz`; reviewable workflow code is under `handoff/r_reference_campaign_20260917/`.
`resources/slurm_job_ledger.csv` retains per-allocation/array-task state, elapsed time and the maximum
available step MaxRSS. Active rows remain explicitly non-final. `slurm_accounting_and_failed_logs.tar.gz`
contains raw accounting steps, job provenance and preserved campaign failure logs; missing measurements
or cancelled-before-start logs are explicitly reported. The remote receipt is separately downloaded and
hash-verified before the local `summary/REMOTE_RECEIPT_UPLOADED.json` success marker is written.
No website or replacement-version notebook is generated, and no old result files are deleted.
''')
    copy(layout)
    campaign_rel=str(OUT.relative_to(ROOT))
    index=STAGE/'README.md'
    index.write_text(f'''# DG-scRNA 实验结果目录

本目录继续使用原 notebook；此前 GBM/PTC 内容保留。本轮 R 参考流程共完成 {summary['expected_analysis_units']} 个分析单元。

- [原结果 notebook：新增实验位于文末](notebooks/dgscrna_results.ipynb)
- [本轮中文结果与解释]({campaign_rel}/summary/RESULTS_ZH.md) · [English report]({campaign_rel}/summary/RESULTS_AND_INTERPRETATION.md)
- [完整流程与 ablation 节点：PDF]({campaign_rel}/summary/workflow_decision_tree.pdf) · [可编辑 SVG]({campaign_rel}/summary/workflow_decision_tree.svg)
- [数据集与 marker 组织清单]({campaign_rel}/summary/marker_context_roster.csv)
- [所有最终注释指标与无 DL 对照]({campaign_rel}/summary/all_annotation_metrics.csv.gz) · [聚类指标]({campaign_rel}/summary/all_clustering_metrics.csv)
- [PTC 成对 ablation 比较]({campaign_rel}/summary/PTC_paired_ablation_changes.csv.gz) · [batch 混合与标签保留诊断]({campaign_rel}/verification/batch_biology/INTERPRETATION.md)
- [逐数据集结果与图]({campaign_rel}/benchmark/) · [PTC ablation 结果与图]({campaign_rel}/PTC_ablation/) · [PTC 原始 checkpoint 对照]({campaign_rel}/PTC_archived_CCA2000/)
- [SLURM 运行账本]({campaign_rel}/resources/slurm_job_ledger.csv) · [目录及归档说明]({campaign_rel}/summary/DELIVERY_LAYOUT.md)
- [交付文件 SHA256 清单]({campaign_rel}/summary/DELIVERY_MANIFEST.json) · [完成后的交付核验记录]({campaign_rel}/summary/DELIVERY_RECEIPT.json)

每个分析目录的 `terminal_results_and_audits.tar.gz` 包含全部初始/最终逐细胞标签、置信度、训练记录、聚类、DEG 和审计资料；在对应目录解压即可恢复明细树。模型权重、概率 NPZ 和大矩阵留在 HPC，并提供原路径及校验值索引。

最终 DG-scRNA 指标使用 DL/refinement 端点；结构性无法训练与无需训练的状态单独保留。历史标签一致性、TCR 代理指标与 reviewer 的 curated labels 各自报告，不能混作同一真值。HCL 按 59 个原组织组运行。论文历史权重与 Accuracy 来源仍未恢复，本轮结果不改变这一限制。
''')
    files.append('README.md')
    files=sorted(set(files))
    manifest=[dict(path=rel,bytes=(STAGE/rel).stat().st_size,sha256=sha(STAGE/rel)) for rel in files]
    manifest_path=PARENT/'R_reference_campaign_20260917_delivery_manifest.json'
    manifest_path.write_text(json.dumps(dict(remote=REMOTE,status='staged',n_files=len(files),
       n_upload_files=len(files)+1,total_bytes=sum(r['bytes'] for r in manifest),files=manifest,
       manifest_scope='Payload hashes; the manifest itself and the later receipt are excluded to avoid self-referential hashes.',
       delivered_manifest=str(OUT.relative_to(ROOT)/'summary/DELIVERY_MANIFEST.json'),
       job=os.environ['SLURM_JOB_ID']),indent=2)+'\n')
    delivered_manifest=OUT/'summary/DELIVERY_MANIFEST.json'
    shutil.copy2(manifest_path,delivered_manifest);copy(delivered_manifest)
    listpath=PARENT/'R_reference_campaign_20260917_upload_files.txt';listpath.write_text('\n'.join(sorted(set(files)))+'\n')
    print('STAGED',len(files),'files including delivery manifest',sum(r['bytes'] for r in manifest)+manifest_path.stat().st_size,'bytes',flush=True)

if __name__=='__main__':run()
