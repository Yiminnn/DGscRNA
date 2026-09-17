"""Report baseline diagnosis in the existing notebook and stage the authorized folder."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import shutil
import sys

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
BASE = ROOT / 'results/hvg_ptc_20260916_v1'
OUT = BASE / 'ptc_paper_baseline'
DELIVERY = BASE / 'onedrive_existing_results_20260916'
STAGE = DELIVERY / 'GBM_PTC_results_20260916'
REMOTE = 'onedrive:work_od/share/dgscrna_GSE274546_TKU3186/01_report/GBM_PTC_results_20260916'


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def run():
    assert os.environ.get('SLURM_JOB_ID')
    import pandas as pd
    diagnosis = OUT / 'recovery_diagnosis'
    summary = json.loads((diagnosis / 'initialization_replicates/SUMMARY.json').read_text())
    routes = pd.read_csv(diagnosis / 'initialization_replicates/route_summary.csv')
    metrics = pd.read_csv(diagnosis / 'initialization_replicates/metrics_all_seeds.csv')
    assert len(routes) == 12 and summary['no_seed_selected']
    assert routes.loc[routes.model_seed.eq(42), 'seed42_control_weights_probabilities_terminal_exact'].eq(True).all()
    seeds = routes[routes.model_seed.ne(42)]
    overall = metrics[metrics.scope.eq('Overall') & metrics.model_seed.ne(42)]
    ranges = summary['ranges']
    audit = json.loads((diagnosis / 'manifest.json').read_text())
    assert audit['source_comparisons'] == 0
    text = '''# PTC 为什么不一致，怎样恢复

原始发表用的最终注释已经找回，并与 Sup 逐细胞核对一致。尚未恢复的是重新训练得到旧模型全部逐细胞输出，以及 V16 Accuracy 行的计算来源。这三件事需要分别验收，不能笼统说“原始结果丢失，所以只能等权重”。

## 已有证据能解释到哪里

1. **此前有过流程选错，已经纠正。** 原始 notebook 的 NMT 是 Thyroid / Seurat / none / 最终 DL，TTU 是 Pubmed34663816 / UMAP-HDBSCAN / mean / 最终 DL；原始 CCA 将八个样本一起整合，再按两组选择不同输出。新建两个独立 CCA 对象会改变特征集，曾使 NMT 丢失 CD3D。当前核验已经改用原始 checkpoint，不能再把这一旧错误当作剩余 DL 差异的已证实原因。
2. **原程序没有固定模型初始化。** `source.py` 的 `Generator().manual_seed(42)` 只用于训练/验证集划分。模型随后随机初始化，原函数顺序训练多个注释列，已有 CSV 还会触发跳过。这意味着在没有全局随机状态记录时，列的执行顺序、断点继续情况都可能改变初始化；这是由源码支持的机制推断，不是历史运行状态已经找回的证明。当前显式使用模型 seed42，不能据此认定复现了旧初始化。
3. **当前实现差异已被限定排除。** 原始 loom 与重训矩阵的全部值、细胞和基因顺序一致；全2000基因 R 检验与计算优化所得结果一致；在相同固定输入、seed42和当前环境下，直接执行原始完整 DL 函数与重建实现的权重、概率、终端标签全一致。这排除了当前 wrapper 或该计算优化在此条件下造成差异，未证明历史初始训练标签、软件环境和随机状态也相同。
4. **剩余差异集中在 DL 补注释。** seed42 的 NMT 有16/48,255、TTU有1,641/44,149个终端标签不同。TTU其中1,315个旧版已注释细胞变为 Unknown，912个原为B细胞。这1,315个细胞的置信度中位数为0.685858，696个低于0.7，因此不能解释为全都恰好落在0.9阈值附近；降低阈值不能恢复旧模型。
5. **仍有一个重要输入缺口。** 已检查原始 cohort 的 CSV、loom 和其他归档注释文件，未直接找回所选两条路线的历史 pre-DL 初始标签列。目前“初始已知标签与旧最终标签不冲突”，不足以证明旧训练细胞集合和标签完全一致。其他样本数的注释文件没有按位置强行拼接。

## 仅改变初始化的 SLURM 诊断

预先固定模型 seeds=0、1、2、3、4，两条路线各跑一次；seed42单独作为实现一致性对照。所有表达值、细胞顺序、重建的初始标签、split seed42、10epochs、训练设置和0.90阈值保持相同。12次拟合全部完成，seed42的权重、概率和最终标签再次完全对上已有重训。没有挑选最佳 seed，也没有用 Sup 最终标签重新训练。

| 模型 seed | NMT 不同细胞 | TTU 不同细胞 | 合并 F1（原稿 class0 定义） | 合并二值 AUC |
|---|---:|---:|---:|---:|
'''
    for seed in [0, 1, 2, 3, 4, 42]:
        rr = routes[routes.model_seed.eq(seed)].set_index('group')
        mm = metrics[metrics.model_seed.eq(seed) & metrics.scope.eq('Overall')].iloc[0]
        label = '42（复核）' if seed == 42 else str(seed)
        text += f"| {label} | {int(rr.loc['NMT','n_mismatches'])} | {int(rr.loc['TTU','n_mismatches'])} | {mm.F1_class0:.6f} | {mm.AUC_binary:.6f} |\n"
    text += ('\n五个预设种子的观察范围为 F1 '
             f"{ranges['F1_class0']['min']:.6f}–{ranges['F1_class0']['max']:.6f}，AUC "
             f"{ranges['AUC_binary']['min']:.6f}–{ranges['AUC_binary']['max']:.6f}。"
             '归档原输出为 F1=0.9518588、AUC=0.9391770。这些范围不是置信区间，不能据此认定全部历史差异仅由初始化导致。\n\n')
    text += '| 分组 | 五种子始终一致的细胞 | 五种子均未恢复原标签的细胞 |\n|---|---:|---:|\n'
    for group, counts in summary['per_group_stability'].items():
        text += f"| {group} | {counts['n_all_5_seeds_agree']}/{counts['n_cells']} | {counts['n_no_seed_matches_historical']} |\n"
    text += '''
## Accuracy 为什么是另一个问题

原标签配原始验证向量已恢复全部8项DG F1/AUC，原稿 F1 的正类为0，且原 T 类映射包含 NK。用完全相同的预测与真值重算普通二分类 Accuracy 是0.941182，V16写0.8922；同一二分类向量下 F1=0.9519 与这个 Accuracy 数学上不相容。其他评价口径、分母或录入错误都有待来源核实，不能据此断言是哪一种。找回的早期 `tcr/scripts/paper.md` 没有 Accuracy 行，不能提供 V16 的来源。

## 怎样恢复

1. **恢复论文原结果：已经可做且主要对账已完成。** 固定 `annotations.xlsx`（S2）和 `data_with_validation.csv`（S3）及其原始 label 简化规则，按 sample+barcode 连接。原始最终注释和8项DG F1/AUC已经恢复；重建论文图表应明确标注使用归档最终输出。保留两个 Sup 终点各自的含义。
2. **恢复旧模型的精确推理：优先需要原始权重，或足够恢复旧训练过程的完整状态。** 最终标签不能唯一反推出权重。剩余取证重点是旧 pre-DL 标签/训练池、当时实际执行的 DL 代码与调用、模型随机状态及环境，而不是继续修改最终标签或把seed42当作旧种子。没有这些记录，不能保证重训到92,404个逐细胞全部相同。
3. **建立以后能重复运行的基线：当前输入、初始标签、split、权重和概率已保存。** 在历史输入缺口核实后，以固定输入和明确的多种子结果对照原论文指标，分别报告逐细胞一致性、各类 Unknown、TCR 指标及运行变异。复现验收不能临时改成“挑一个最像旧结果的seed”；下一步消融要用同一组预设种子成对比较。当前诊断不自动视为通过原始基线门槛，PTC消融继续暂停。
4. **恢复 Accuracy 的可解释性：查清V16定义与来源；若最终确认采用同一二分类普通accuracy，应有依据地更正为可重算的0.941182并同步定义。** 如果它来自另一个终点，应明确该终点、分母和代码。当前未修改原稿数值，也不调模型去追0.8922。

原始输出、现在可重跑的模型、初始化诊断分别留存。原 notebook 更新原有PTC说明单元，所有GBM分析和图保持完整；新增诊断同步到已授权 OneDrive 原目录。
'''
    report = OUT / 'RESTORATION_DIAGNOSIS_ZH.md'
    report.write_text(text)
    status = OUT / 'STATUS.md'
    marker = '\n\n## Additional restoration diagnosis (2026-09-17)\n'
    s = status.read_text().split(marker)[0]
    s += marker + ('A prespecified five-model-seed diagnostic (0–4, plus seed42 parity controls) completed for both routes. '
        'See [restoration diagnosis](RESTORATION_DIAGNOSIS_ZH.md) and `recovery_diagnosis/initialization_replicates/SUMMARY.json`. '
        'All twelve fits and seed42 parity controls passed. Historical pre-DL training-label identity is still not established; '
        'this diagnostic does not establish a unique historical cause or automatically pass the original-workflow gate.\n')
    status.write_text(s)
    state_path = OUT / 'RECONSTRUCTION_STATUS.json'
    state = json.loads(state_path.read_text())
    state['restoration_diagnosis'] = dict(report=str(report), summary=summary, historical_initial_labels_directly_recovered=False)
    state_path.write_text(json.dumps(state, indent=2) + '\n')
    import correct_notebook
    correct_notebook.run()
    task = ROOT / 'handoff/ptc_paper_baseline_20260916/CURRENT_TASK.md'
    task_text = task.read_text().split(marker)[0]
    task_text += marker + ('The user asked why reruns differ and how to restore. Fixed-input initialization diagnostics '
        '(array7340109, summary7340110; five seeds plus seed42 controls per route) completed. '
        'No seed was selected. Read `RESTORATION_DIAGNOSIS_ZH.md` and the new diagnosis manifests before attributing the remaining gap to randomness. '
        'Exact original final labels are already restored; historical initial-training-label equality and V16 Accuracy provenance remain unresolved. '
        'Check the separate `PTC_recovery_diagnosis_upload_receipt.json` for this addition’s delivery.\n')
    task.write_text(task_text)
    notebook_manifest = OUT / 'notebook_correction_manifest.json'
    delivery_state_path = BASE / 'EXPERIMENT_DELIVERY.json'
    delivery_state = json.loads(delivery_state_path.read_text())
    notebook_info = json.loads(notebook_manifest.read_text())
    delivery_state['notebook_sha256'] = notebook_info['current_sha256']
    delivery_state['notebook_updated_at'] = notebook_info['updated_utc']
    delivery_state['PTC_restoration_diagnosis'] = str(report)
    delivery_state_path.write_text(json.dumps(delivery_state, indent=2) + '\n')
    paths = [report, status, state_path, task, notebook_manifest, delivery_state_path, ROOT/'notebooks/dgscrna_results.ipynb']
    paths.extend(p for p in diagnosis.rglob('*') if p.is_file())
    script_dir = ROOT / 'handoff/ptc_paper_baseline_20260916'
    paths.extend(script_dir/name for name in ['diagnose_recovery_inputs.py','diagnose_recovery_inputs.sbatch',
        'diagnose_initialization.py','diagnose_initialization.sbatch','summarize_initialization.sbatch',
        'correct_notebook.py','finalize_recovery_diagnosis.py','finalize_recovery_diagnosis.sbatch'])
    records = []
    for source in sorted(set(paths)):
        relative = source.relative_to(ROOT)
        target = STAGE / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        checksum = sha(source)
        assert checksum == sha(target)
        records.append(dict(path=str(relative), bytes=source.stat().st_size, sha256=checksum))
    (DELIVERY / 'PTC_recovery_diagnosis_upload_files.txt').write_text('\n'.join(r['path'] for r in records) + '\n')
    staging = dict(job=os.environ['SLURM_JOB_ID'], staged_utc=datetime.now(timezone.utc).isoformat(),
                   remote=REMOTE, files=records, n_files=len(records), total_bytes=sum(r['bytes'] for r in records))
    (DELIVERY / 'PTC_recovery_diagnosis_staging.json').write_text(json.dumps(staging, indent=2) + '\n')
    print(json.dumps({k:v for k,v in staging.items() if k!='files'}, indent=2), flush=True)


def receipt():
    staging = json.loads((DELIVERY / 'PTC_recovery_diagnosis_staging.json').read_text())
    log = DELIVERY / 'PTC_recovery_diagnosis_remote_check.log'
    assert '0 differences found' in log.read_text()
    staging['verified_utc'] = datetime.now(timezone.utc).isoformat()
    staging['verification'] = 'rclone check --download --one-way completed successfully; zero differences'
    staging['verification_log_sha256'] = sha(log)
    (DELIVERY / 'PTC_recovery_diagnosis_upload_receipt.json').write_text(json.dumps(staging, indent=2) + '\n')
    print('Verified diagnosis delivery:', staging['n_files'], 'files', flush=True)


if __name__ == '__main__':
    receipt() if len(sys.argv)>1 and sys.argv[1]=='receipt' else run()
