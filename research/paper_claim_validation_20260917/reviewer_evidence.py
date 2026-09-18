"""Map the actual PI/reviewer requests to verified artifacts and remaining limits."""
import json
import os
import re
from common import ROOT,OUT,OLD,require_slurm,checked,write_json,sha,utc

def notebook_text(phase):
    dest=OUT/('GBM_full_summary' if phase=='GBM' else 'PTC_summary')
    text=(dest/'REVIEWER_PI_EVIDENCE.md').read_text()
    return re.sub(r'\]\(([^)]+)\)',lambda m:']('+os.path.relpath((dest/m[1]).resolve(),ROOT/'notebooks')+')',text)

def run(phase='GBM'):
    require_slurm()
    assert phase in ['GBM','PTC']
    if phase=='PTC':
        from ptc_followup_common import require_ptc
        require_ptc()
    dest=OUT/('GBM_full_summary' if phase=='GBM' else 'PTC_summary');dest.mkdir(exist_ok=True)
    def link(path,label):return f'[{label}]({os.path.relpath(path,dest)})'
    def status(folder):return '已核验' if checked(OUT/folder) else '待完成'
    core='已核验' if checked(OUT/'summary','aggregate_manifest.json','AGGREGATE_COMPLETE') else '待完成'
    previous=json.loads((OLD/'summary/campaign_summary.json').read_text())
    assert previous['status']=='complete' and previous['evaluated_units']==100
    rows=[
      ('PI08-13：固定 marker 的降维/聚类与逐样本图；R2-2模块消融',core,
       link(OUT/'summary/primary_fixed_marker_24.csv','六档HVG×四路线')+'；'+link(OUT/'summary/workflow_decision_tree.pdf','英文主流程图'),
       '完整原R终端注释；HVG作用按RNA与CCA分支解释；每种聚类有图。最高点只适用于列明的搜索范围。'),
      ('PI08-04/08-13：不同组织及外部研究 marker',status('marker_evidence_summary'),
       link(OUT/'marker_evidence_summary/library_source_assay_coverage_audit.csv','来源/技术/覆盖审计'),
       '16库完整保留；参与原作者标签构建的CARE/BrainAtlas及包含它们的union不参与主选择。CellMarker组织子库不自动等于独立发表reference。'),
      ('R2-1、R3-B1：公共多类数据、逐类评价','既有原R campaign已完成',
       link(OLD/'summary/cohort_sizes.csv','11个reviewer数据集')+'；'+link(OLD/'summary/RESULTS_AND_INTERPRETATION.md','原R结果与解释'),
       '100分析单元/12048终端条件；HCL按59组织运行不等于单次60万细胞拟合。不同可支持标签层级分别报告。'),
      ('R3-B2：避免测试标签选marker，给对手可比机会',status('comparison_summary'),
       link(OUT/'comparison_summary/training_patient_choices.csv','训练患者选择')+'；'+link(OUT/'comparison_summary/equal24_DG_scType_comparison.csv','等24配置比较'),
       '固定聚类下13个主marker库，另报DG/scType同等24配置选择；SingleR及GNN使用不同reference信息。回顾性患者留出，未声称未见外部测试。'),
      ('R3-B4：全部基线的绝对/相对差及统计',status('comparison_summary'),
       link(OUT/'comparison_summary/paired_patient_comparisons.csv','GBM配对差与增益')+'；'+link(OUT/'PTC_comparator_replay/paired_patient_gains.csv','PTC配对增益'),
       'GBM按患者配对及Holm；PTC保留4患者全部差异，不能由细胞数或种子扩充生物重复。PTC链接完成状态由其manifest决定。'),
      ('R2-3：GNN对手',status('comparison_summary'),
       link(OUT/'comparison_summary/patient_heldout_summary.csv','scDeepSort与其他方法'),
       '已发表Brain预训练GNN；缺乏Malignant输出类的限制保留在分母。不是给DG的MLP虚构图输入。'),
      ('R2-5及R3超参数/学习曲线',status('controls_summary'),
       link(OUT/'controls_summary/MLP_all_metrics.csv','MLP控制')+'；'+link(OUT/'controls_summary/representation_all_metrics.csv','表示/聚类控制'),
       '原参数、种子、维数、minPts/resolution、MLP宽度/epoch；0.7/0.9共用保存概率。不存在于源码的alpha/epsilon通过方法文字纠正。参数敏感性限3个事先选定pilot。'),
      ('R1-2、R2-4、R3-C3：Unknown表达和错误传播',status('unknown_expression'),
       link(OUT/'unknown_expression/all_patient_paired_gene_contrasts.csv.gz','患者配对表达差')+'；'+link(OUT/'unknown_summary/retained_seed_and_DL_new_errors.csv','种子保留/新增错误'),
       '同样本同原细胞大类内Unknown比较；支持数、效应、CI及全局BH均保留。算法由表达定义，关联不能证明双细胞或新群体。'),
      ('R2-6：技术batch与组织生物背景','既有比较已完成',
       link(OLD/'summary/PTC_batch_biology_compact.csv','CCA/NONE/Harmony与生物保留'),
       '比较实际校正流程；不是单纯校正先后顺序效应。原all8与用户指定NMT/TTU分开；共享细胞类型内考察混合。'),
      ('R2-9、R3-C1：>100k及多次耗时/内存',status('scalability_summary'),
       link(OUT/'scalability_summary/resource_repeated_mean_SD.csv','3方法×5规模×3运行')+'；'+link(OUT/'scalability_summary/RESOURCE_INTERPRETATION.md','资源口径'),
       '单次最大120k，独立进程及隔离模型缓存；实际CPU分配、节点、wall time与job/step MaxRSS。CPU运行的GPU显存不适用；OS缓存不受控。'),
      ('R3-B3：核验SignacX',status('PTC_comparator_replay'),
       link(OUT/'PTC_comparator_replay/ENDPOINT_CORRECTIONS.md','缓存逐细胞重评与解释修正'),
       'CellStates显式T与CellTypes TNK分别评价。原S3有标签但无淋巴类；不再声称运行崩溃，也不把未隔离的版本/graph因素当确定原因。'),
      ('PTC复现锚点、marker×DL、患者选择、种子及marker保留',status('PTC_summary'),
       link(OUT/'PTC_summary/claim_summary_by_group.csv','分组实际结论')+'；'+link(OUT/'PTC_summary/paired_patient_contrasts.csv','配对机制对照'),
       '保留NMT-SNN、TTU-UMAP/HDBSCAN及不同原marker。严格T/生产性TCR、历史non-T F1/AUC与保存标签一致率分开。'),
      ('R3-A：模型framing与方法源码一致','提供方法修订文字',
       link(OUT/'GBM_full_summary/METHODS_CORRECTIONS.md','源码支持的方法文字'),
       '现有模型为表达MLP refinement；图只参与上游聚类。稿件本体未被改写，不宣称这些文字已经投稿。'),
      ('R3-C2、R2-7：新亚群、转化/治疗结论','需要限定主张',
       link(ROOT/'paper/comments.md','reviewer原文'),
       '本轮计算不能代替独立功能验证或扩大4患者生物学样本。marker/QC/TCR关联不能证明新细胞群或治疗靶点。'),
      ('GEO/Zenodo、湿实验参数及原Accuracy来源','未由本轮计算解决',
       link(ROOT/'paper/comments.md','原始要求')+'；'+link(ROOT/'handoff/paper_claim_validation_20260917/README.md','运行代码范围'),
       '代码仅更新指定branch；未创建永久DOI或操作未经确认账户。原始Accuracy行计算来源及历史权重仍缺失，近似复现范围明确。')]
    text='''# PI / reviewer evidence index

下表将实际要求对应到结果，不把原始方法的有效性、内部最高配置和方法间优越性混为一项。
本轮沿用原notebook与原OneDrive目录；没有创建网页。原数据、旧结果和反例保留。
状态依据完成manifest生成；待完成链接不是已经通过的实验。

| 要求 | 当前证据状态 | 结果入口 | 可支持的解释与限制 |
|---|---|---|---|
'''
    text+='\n'.join('| '+' | '.join(row)+' |' for row in rows)+'\n'
    text+='\nPI来源：'+link(ROOT/'handoff/meeting_0813/transcript_full.txt','08-13转录')+'；'+link(ROOT/'handoff/meeting_0804/ACTION_ITEMS_0804.md','08-04任务')+'；'+link(ROOT/'handoff/meeting_0730/ACTION_ITEMS_0730.md','07-30任务')+'。\n'
    (dest/'REVIEWER_PI_EVIDENCE.md').write_text(text)
    write_json(dest/'reviewer_evidence_manifest.json',dict(phase=phase,requirements=len(rows),
        status='evidence_index_with_remaining_limits',old_campaign_sha256=sha(OLD/'summary/campaign_summary.json'),
        report_sha256=sha(dest/'REVIEWER_PI_EVIDENCE.md'),job=os.environ['SLURM_JOB_ID'],completed_at=utc()))

if __name__=='__main__':run()
