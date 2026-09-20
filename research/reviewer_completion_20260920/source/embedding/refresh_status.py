"""Read only A1 completion manifests, never matrices or scientific metrics."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
CAMPAIGN = ROOT / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
OUT = CAMPAIGN / 'embedding'


def checked(directory, manifest='manifest.json', flag='COMPLETE'):
    try:
        return (directory / flag).read_text().strip() == hashlib.sha256((directory / manifest).read_bytes()).hexdigest()
    except OSError:
        return False


protocol = json.loads((CAMPAIGN / 'protocol/embedding.json').read_text())
geometry, fitted, evaluated, in_progress = [], [], [], []
for sample in protocol['samples']:
    for budget in protocol['budgets']:
        base = OUT / sample / budget
        if checked(base / 'geometry'):
            geometry.append(f'{sample}/{budget}')
        for space in protocol['spaces']:
            directory = base / space
            if checked(directory, 'fit_manifest.json', 'FIT_COMPLETE'):
                fitted.append(f'{sample}/{budget}/{space}')
                if checked(directory / 'evaluation') and checked(directory / 'figures'):
                    evaluated.append(f'{sample}/{budget}/{space}')
            elif (directory / 'config.json').exists():
                try:
                    config = json.loads((directory / 'config.json').read_text())
                    conditions = config.get('conditions', [])
                    if conditions:
                        scored = sum(checked(Path(c['dest']), 'score_manifest.json', 'SCORE_COMPLETE') for c in conditions)
                        terminal = sum(checked(Path(c['dest']) / 'terminal/L00_mean', 'terminal_manifest.json', 'TERMINAL_COMPLETE') for c in conditions)
                        in_progress.append(f'{sample}/{budget}/{space}：评分 {scored}/13；终端 DL {terminal}/13')
                except (ValueError, OSError):
                    pass
jobs = json.loads((OUT / 'jobs.json').read_text()) if (OUT / 'jobs.json').exists() else []
status = dict(stage='A1_GBM', work_package='A', status='running',
    updated_at=datetime.now(timezone.utc).isoformat(), jobs=jobs,
    summary='七种表示与三类聚类；计数仅含评分、终端DL、评价和全部图均有完成凭据的候选。',
    completed=len(evaluated) * 13, remaining=22022 - len(evaluated) * 13,
    details=[f'原R几何输入导出：{len(geometry)}/242 个 sample×HVG 单元',
             f'终端拟合完成：{len(fitted)}/1694 个表示单元',
             f'评价和图完成：{len(evaluated)}/1694 个表示单元',
             '每个表示单元包含6个K-means、6个GMM及1个HDBSCAN候选；全队列选K及科学验收另记。',
             *in_progress],
    evidence=[str(p.relative_to(ROOT)) for p in [CAMPAIGN / 'protocol/embedding.json', CAMPAIGN / 'protocol/embedding_selection.json']],
    whole_work_package_A_complete=False)
path = OUT / 'status.json';tmp = path.with_name('status.json.' + str(os.getpid()) + '.tmp')
tmp.write_text(json.dumps(status, ensure_ascii=False, indent=2) + '\n');tmp.replace(path)
print(json.dumps(dict(geometry=len(geometry), fitted=len(fitted), evaluated=len(evaluated))))
