#!/usr/bin/env python3
"""Overlay live execution facts on the existing compact reviewer page.

This never changes review_data.json, existing verdicts, or source provenance.
Run via SLURM. --watch 45 updates the local HTML snapshot until --stop-file
exists or progress.json has watcher_stop=true; it does not submit/cancel jobs.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import html
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
CAMPAIGN = ROOT / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
PAGE = ROOT / 'results/hvg_ptc_20260916_v1/paper_review_text_20260920/index.html'
WEB = ROOT / 'handoff/reviewer_completion_20260920/web'
START = '<!-- REVIEWER_EXECUTION_PROGRESS_BEGIN -->'
END = '<!-- REVIEWER_EXECUTION_PROGRESS_END -->'
PACKAGE_TITLES = {
    'A': 'GBM 流程', 'B': '参数与训练', 'C': '公平方法比较',
    'D': 'Marker 与生物证据', 'E': 'Batch 对照', 'F': '一致性与资源',
    'G': '稿件与回复', 'H': '来源与复现',
}
STATUS_LABELS = {
    'not_started': '未开始', 'planned': '待执行', 'pending': '待执行',
    'running': '进行中', 'in_progress': '进行中', 'preparing': '准备中',
    'implementing': '实现中', 'parity_submitted': '锚点核验已提交',
    'partial': '部分完成', 'blocked': '受阻', 'failed': '失败·待修复',
    'needs_retry': '待重跑', 'awaiting_validation': '待验收',
    'completed': '完成·待验收', 'complete': '完成·待验收',
    'validated': '已验收', 'validated_complete': '已验收',
    'reused': '复用已有证据', 'awaiting_source': '待来源',
    'awaiting_publication': '待发布', 'stopped': '已停止',
}
JOB_LABELS = {
    'RUNNING': '运行', 'PENDING': '排队', 'CONFIGURING': '配置中',
    'COMPLETING': '收尾', 'COMPLETED': '作业结束·待验收',
    'FAILED': '失败', 'OUT_OF_MEMORY': 'OOM·待重跑', 'TIMEOUT': '超时',
    'CANCELLED': '取消', 'NODE_FAIL': '节点故障', 'PREEMPTED': '抢占',
}


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def atomic_write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f'.{os.getpid()}.tmp')
    tmp.write_text(text)
    tmp.replace(path)


def read_json(path, diagnostics):
    try:
        value = json.loads(path.read_text())
        if not isinstance(value, dict):
            raise ValueError('expected JSON object')
        return value
    except FileNotFoundError:
        diagnostics.append(f'尚无状态文件：{path.relative_to(ROOT)}')
    except (ValueError, OSError) as exc:
        diagnostics.append(f'状态暂不可读：{path.relative_to(ROOT)} ({type(exc).__name__})')
    return {}


def as_list(value):
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def plain(value):
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False)
    return '' if value is None else str(value)


def esc(value):
    return html.escape(plain(value), quote=True)


def packages_from(document):
    source = document.get('work_packages', document.get('packages', []))
    if isinstance(source, dict):
        return [{'id': key, **(value if isinstance(value, dict) else {'summary': value})}
                for key, value in source.items()]
    return [p for p in as_list(source) if isinstance(p, dict)]


def job_ids(package):
    ids = []
    for job in as_list(package.get('jobs', package.get('slurm_jobs', []))):
        job_id = plain(job.get('job_id', job.get('id', ''))) if isinstance(job, dict) else plain(job)
        if re.fullmatch(r'\d+(?:_[\d,\[\]%\-]+)?(?:\.\d+)?', job_id):
            ids.append(job_id)
    for key in ('job_id', 'slurm_job_id'):
        job_id = plain(package.get(key, ''))
        if re.fullmatch(r'\d+(?:_[\d,\[\]%\-]+)?(?:\.\d+)?', job_id):
            ids.append(job_id)
    return sorted(set(ids))


def scheduler_snapshot(ids):
    """Only read the campaign's explicit jobs, including completed array tasks."""
    if not ids:
        return {'checked_at': datetime.now(timezone.utc).isoformat(), 'jobs': {}, 'errors': []}
    base_ids = sorted({re.split(r'[_\.]', i)[0] for i in ids})
    jobs, errors = {}, []
    for command, source in [
        (['sacct', '--array', '-n', '-X', '-P', '-j', ','.join(base_ids),
          '--format=JobID%80,State,ExitCode,Elapsed'], 'sacct'),
        (['squeue', '--array', '-h', '-j', ','.join(base_ids), '-o', '%i|%T|%M|%R'], 'squeue'),
    ]:
        try:
            result = subprocess.run(command, capture_output=True, text=True, timeout=18)
            if result.returncode:
                errors.append(f'{source} 查询失败；保留未知状态')
                continue
            for line in result.stdout.splitlines():
                parts = line.strip().split('|')
                if len(parts) < 4:
                    continue
                jid, state = parts[0], parts[1].split()[0].rstrip('+')
                jobs[jid] = {'state': state, 'source': source,
                             ('exit_code' if source == 'sacct' else 'elapsed'): parts[2],
                             ('elapsed' if source == 'sacct' else 'reason_or_node'): parts[3]}
        except (OSError, subprocess.TimeoutExpired):
            errors.append(f'{source} 查询超时或不可用；未判断为停止')
    return {'checked_at': datetime.now(timezone.utc).isoformat(), 'jobs': jobs, 'errors': errors}


def states_for(job_id, snapshot):
    jobs = snapshot['jobs']
    found = [record['state'] for jid, record in jobs.items()
             if jid == job_id or jid.startswith(job_id + '_')]
    return Counter(found) if found else Counter({'UNKNOWN': 1})


def job_text(job_id, snapshot):
    counts = states_for(job_id, snapshot)
    return job_id + ' ' + ' / '.join(
        JOB_LABELS.get(state, '状态未知' if state == 'UNKNOWN' else state)
        + (f' ×{count}' if count > 1 else '') for state, count in sorted(counts.items()))


def evidence_items(package):
    return as_list(package.get('evidence', package.get('artifacts', [])))


def state_label(package):
    raw = plain(package.get('status', 'not_started'))
    value = raw.lower()
    if package.get('verified') is True and package.get('verification') and evidence_items(package):
        path = (ROOT / package['verification']).resolve()
        if path.is_relative_to(ROOT):
            try:
                if json.loads(path.read_text()).get('status') == 'passed':
                    return '已验收', 'validated'
            except (OSError, ValueError):
                pass
    # Do not promote a process completion or undocumented claim to validated.
    if value in ('validated', 'validated_complete'):
        if not package.get('validation') or not evidence_items(package):
            return '待验收（验收证据未登记）', 'awaiting_validation'
    return STATUS_LABELS.get(value, raw), value


def evidence_html(item):
    raw = plain(item.get('path', item.get('url', ''))) if isinstance(item, dict) else plain(item)
    label = plain(item.get('label', raw)) if isinstance(item, dict) else raw
    if not raw:
        return ''
    candidate = (ROOT / raw).resolve() if not raw.startswith('/') else Path(raw).resolve()
    if candidate.is_relative_to(ROOT) and candidate.exists():
        href = os.path.relpath(candidate, PAGE.parent)
        return f'<a href="{esc(href)}">{esc(label)}</a>'
    return f'<span>{esc(label)}（路径待核验）</span>'


def detail_html(package, snapshot, req_ids):
    rows = []
    for key, title in [('summary', '当前'), ('completed', '已处理'), ('remaining', '待处理'), ('details', '详情'), ('next', '下一步'),
                       ('blockers', '限制'), ('scope', '范围'), ('validation', '验收'), ('verification', '验收文件'), ('updated_at', '登记时间')]:
        if package.get(key):
            value = '\n'.join(map(plain, as_list(package[key])))
            rows.append(f'<dt>{title}</dt><dd>{esc(value)}</dd>')
    ids = job_ids(package)
    if ids:
        rows.append('<dt>SLURM</dt><dd>' + '<br>'.join(esc(job_text(j, snapshot)) for j in ids) + '</dd>')
    requirements = [r for r in as_list(package.get('requirements', [])) if r in req_ids]
    if requirements:
        links = ' · '.join(f'<button class="link" data-action="select-requirement" data-requirement="{esc(r)}">{esc(r)}</button>' for r in requirements)
        rows.append(f'<dt>对应要求</dt><dd>{links}</dd>')
    evidence = [evidence_html(x) for x in evidence_items(package)]
    if any(evidence):
        rows.append('<dt>已登记材料</dt><dd>' + '<br>'.join(filter(None, evidence)) + '</dd>')
    return '<dl>' + ''.join(rows) + '</dl>' if rows else '<p>尚无新增执行记录；已有科学证据见下方要求列表。</p>'


def worker_units(worker_name, worker):
    units = packages_from(worker)
    if units:
        return units
    target = worker.get('work_package', worker.get('package'))
    if not target and re.match(r'^[A-H](?:_|$)', plain(worker.get('stage'))):
        targets = []
        for token in worker['stage'].split('_'):
            if token not in PACKAGE_TITLES:
                break
            targets.append(token)
        return [{**worker, 'id': key} for key in targets]
    if not target and worker_name == 'controls' and worker:
        target = 'B'
    return [{**worker, 'id': target}] if target else []


def render_block(progress, workers, snapshot, diagnostics, review):
    req_ids = {r['id'] for r in review['requirements']}
    source = {plain(p.get('id', '')).upper(): p for p in packages_from(progress)}
    rows = []
    for key, default_title in PACKAGE_TITLES.items():
        package = source.get(key, {'id': key, 'title': default_title, 'status': 'not_started'})
        label, state = state_label(package)
        combined_job_ids = set(job_ids(package))
        children = []
        for worker_name, worker in workers:
            for unit in worker_units(worker_name, worker):
                target = plain(unit.get('id', unit.get('work_package', unit.get('package', '')))).upper()
                if target != key:
                    continue
                combined_job_ids.update(job_ids(unit))
                child_label, _ = state_label(unit)
                scope = plain(unit.get('stage', worker.get('stage', worker_name)))
                children.append(f'<p class="execution-worker">子任务 {esc(scope)}：{esc(child_label)}（不代表整个工作包验收）</p>' + detail_html(unit, snapshot, req_ids))
        brief_jobs = ' · '.join(job_text(j, snapshot) for j in sorted(combined_job_ids))
        rows.append(f'<details class="execution-row" id="execution-{key}" data-package="{key}" data-status="{esc(state)}">'
                    f'<summary><span class="execution-id">{key}</span> {esc(package.get("title", default_title))} '
                    f'<span class="execution-state">[{esc(label)}]</span>'
                    + (f'<span class="execution-job"> · {esc(brief_jobs)}</span>' if brief_jobs else '')
                    + '</summary><div class="execution-detail">' + detail_html(package, snapshot, req_ids)
                    + ''.join(children) + '</div></details>')
    now = datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S %Z')
    notes = [*as_list(progress.get('notes')), *diagnostics, *snapshot['errors']]
    notes_html = ''.join(f'<p>{esc(note)}</p>' for note in notes)
    return START + '''
<style id="reviewer-execution-style">
.execution-progress{padding:5px 0 9px;font-size:12px}.execution-heading{display:flex;align-items:baseline;gap:10px;flex-wrap:wrap}.execution-heading h2{font-size:13px;margin:0}.execution-stamp,.execution-hint{font-size:10px;color:#777}.execution-hint{margin:1px 0 3px}.execution-row{border:0;margin:0;scroll-margin-top:60px}.execution-row>summary{cursor:pointer;list-style:none;position:relative;min-height:24px;line-height:24px;padding-left:13px}.execution-row>summary::-webkit-details-marker{display:none}.execution-row>summary:before{content:'›';position:absolute;left:0;color:#777}.execution-row[open]>summary:before{transform:rotate(90deg)}.execution-id{color:#777;font-size:11px;margin-right:4px}.execution-state{font-size:11px;color:#526b77;margin-left:5px}.execution-row[data-status="failed"] .execution-state,.execution-row[data-status="blocked"] .execution-state{color:#ac4a39}.execution-row[data-status="validated"] .execution-state,.execution-row[data-status="validated_complete"] .execution-state{color:#28704a}.execution-job{font-size:10px;color:#777;overflow-wrap:anywhere}.execution-detail{padding:2px 0 8px 24px;font-size:12px;max-width:1000px}.execution-detail dl{display:grid;grid-template-columns:64px minmax(0,1fr);gap:3px 10px;margin:3px 0}.execution-detail dt{color:#777}.execution-detail dd{margin:0;white-space:pre-wrap;overflow-wrap:anywhere}.execution-detail p{margin:3px 0}.execution-worker{font-weight:550}.execution-notes{font-size:10px;color:#777;margin-top:3px}.execution-notes summary{cursor:pointer}.execution-notes p{margin:4px 0;overflow-wrap:anywhere}.execution-heading button{font-size:10px;color:#52734e;text-decoration:underline}.execution-links{font-size:10px;margin:1px 0 3px;overflow-wrap:anywhere}@media(max-width:760px){.execution-detail{padding-left:12px}.execution-detail dl{grid-template-columns:54px minmax(0,1fr);font-size:11px}.execution-job{font-size:9px}}@media print{.execution-progress{display:none!important}}
</style>
<section class="execution-progress" id="reviewer-execution" aria-label="Reviewer 补齐执行进度">
<div class="execution-heading"><h2>执行进度</h2><time class="execution-stamp">''' + esc(now) + '''</time></div>
<p class="execution-hint">点行展开 · 每45秒更新文件；Paseo中关闭后重新打开查看 · 作业结束与科学验收分开记录</p>
<p class="execution-links"><a href="../../../handoff/reviewer_completion_plan_20260920/PLAN.md">完整计划</a> · <a href="../../../handoff/reviewer_completion_plan_20260920/REVIEWER_COVERAGE.md">35 项验收表</a></p>
''' + ''.join(rows) + '<details class="execution-notes"><summary>更新说明 / 状态来源</summary>' + notes_html + '<p>本区只汇报执行状态；下方 51 项要求、24 项消融、正反证据和出处保持原样。原稿 SignacX 保留 0 T cell，后续核验独立记录。</p></details></section>\n' + END


def update_once():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('HTML rendering must run through SLURM')
    WEB.mkdir(parents=True, exist_ok=True)
    diagnostics = []
    for worker_name in ('embedding', 'no_clustering'):
        refresh = ROOT / 'handoff/reviewer_completion_20260920' / worker_name / 'refresh_status.py'
        registry = CAMPAIGN / worker_name / 'jobs.json'
        if refresh.exists() and registry.exists():
            try:
                result = subprocess.run([sys.executable, str(refresh)], capture_output=True, text=True, timeout=30)
                if result.returncode:
                    diagnostics.append(f'{worker_name} 进度刷新失败，保留上次记录')
            except (OSError, subprocess.TimeoutExpired):
                diagnostics.append(f'{worker_name} 进度刷新超时，保留上次记录')
    progress = read_json(CAMPAIGN / 'progress.json', diagnostics)
    if not progress and (CAMPAIGN / 'progress.json').exists():
        raise RuntimeError('Root progress file temporarily unreadable; keeping previous HTML')
    workers = [(name, read_json(CAMPAIGN / name / 'status.json', diagnostics))
               for name in ('controls', 'evidence', 'manuscript', 'comparison',
                            'comparison/scDeepSort_LogNormalize', 'embedding', 'no_clustering')]
    # These administrative feeds appear once their scoped work has started.
    for name in ('resources', 'embedding/dispatch'):
        if (CAMPAIGN / name / 'status.json').exists():
            workers.append((name, read_json(CAMPAIGN / name / 'status.json', diagnostics)))
    original = PAGE.read_text()
    base = re.sub(re.escape(START) + r'.*?' + re.escape(END), '', original, flags=re.S)
    if not (WEB / 'base_index.html').exists():
        atomic_write(WEB / 'base_index.html', base)
    match = re.search(r'<script id="review-data" type="application/json">(.*?)</script>', base, re.S)
    if not match:
        raise RuntimeError('Existing review metadata missing; refusing replacement')
    review_payload = match.group(1)
    review = json.loads(review_payload)
    assert len(review['requirements']) == 51
    assert sum(len(n['ablations']) for n in review['nodes']) == 24
    ids = set()
    for document in [progress, *(w for _, w in workers)]:
        ids.update(job_ids(document))
        for p in packages_from(document):
            ids.update(job_ids(p))
    snapshot = scheduler_snapshot(sorted(ids))
    block = render_block(progress, workers, snapshot, diagnostics, review)
    if '<main>' not in base:
        raise RuntimeError('Main insertion point missing')
    new_html = base.replace('<main>', '<main>' + block, 1)
    recovered = re.sub(re.escape(START) + r'.*?' + re.escape(END), '', new_html, flags=re.S)
    assert recovered == base
    assert new_html.count('id="reviewer-execution"') == 1
    atomic_write(PAGE, new_html)
    atomic_write(PAGE.with_name('DGscRNA_审阅列表.html'), new_html)
    receipt = {'status': 'rendered', 'rendered_at': datetime.now(timezone.utc).isoformat(),
               'slurm_job_id': os.environ['SLURM_JOB_ID'], 'slurm_step_id': os.environ.get('SLURM_STEP_ID'),
               'page': str(PAGE), 'review_payload_sha256': sha(review_payload),
               'base_html_sha256': sha(base), 'html_sha256': sha(new_html),
               'requirements': 51, 'ablations': 24, 'work_packages': 8,
               'input_status_updated_at': progress.get('updated_at'), 'scheduler': snapshot,
               'diagnostics': diagnostics, 'original_scientific_data_unchanged': True}
    atomic_write(WEB / 'render_receipt.json', json.dumps(receipt, ensure_ascii=False, indent=2) + '\n')
    atomic_write(CAMPAIGN / 'web' / 'render_receipt.json', json.dumps(receipt, ensure_ascii=False, indent=2) + '\n')
    return progress, receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--watch', type=int, default=0, help='Update every N seconds (minimum 15)')
    parser.add_argument('--stop-file', type=Path, default=CAMPAIGN / 'web' / 'STOP')
    args = parser.parse_args()
    while True:
        iteration_started = time.monotonic()
        if args.stop_file.exists():
            break
        try:
            progress, receipt = update_once()
            print(json.dumps({'status': receipt['status'], 'updated_at': receipt['rendered_at'],
                              'jobs': len(receipt['scheduler']['jobs'])}), flush=True)
            if not args.watch or progress.get('watcher_stop'):
                break
        except Exception as exc:
            if not args.watch:
                raise
            print(json.dumps({'status': 'update_failed_previous_html_kept',
                              'error': str(exc)}, ensure_ascii=False), flush=True)
        time.sleep(max(1, max(15, args.watch) - (time.monotonic() - iteration_started)))


if __name__ == '__main__':
    main()
