"""Shared, reference-only postprocessing for accuracy and timing experiments."""
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
from zoneinfo import ZoneInfo

COMPONENTS = ('rho', 'momentum_x', 'momentum_y', 'energy')
RESIDUAL_NOTE = ('HDG Krylov 控制 Schur 残差；恢复后的 DG 残差按 '
                 'max(10×ltol×||R||_(M^-1), allowance) 验收。'
                 '两种方法不能视为每次 Krylov 都采用相同 DG 残差度量。')


def check(condition, message):
    if not condition:
        raise ValueError(message)


def close(a, b):
    return math.isfinite(float(a)) and math.isclose(float(a), float(b), rel_tol=1e-11, abs_tol=1e-14)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Run:
    def __init__(self, root, kind):
        self.root = Path(root).resolve()
        self.inputs = set()
        self.job = self.text(self.root/'job_id.txt').strip()
        check(self.job.isdigit(), 'Expected a numeric job_id.txt')
        self.job_root = self.root/f'job_{self.job}'
        self.data = self.job_root/'data'
        check((self.data/'SUCCESS').is_file(), 'Requested sweep is incomplete')
        self.text(self.data/'SUCCESS')
        self.kind = kind
        self.identity = self.metadata()

    def text(self, path):
        self.inputs.add(path)
        return path.read_text()

    def csv(self, path):
        return list(csv.DictReader(self.text(path).splitlines()))

    def json(self, path):
        return json.loads(self.text(path))

    def metadata(self):
        manifest = self.root/'source.sha256'
        digest = hashlib.sha256(self.text(manifest).encode()).hexdigest()
        base = self.text(self.root/'base_commit.txt').strip()
        exe = self.json(self.data/'executable.json')
        env = self.text(self.job_root/'environment.txt')
        match = re.search(r'^Model name:\s*(.+)$', env, re.MULTILINE)
        result = dict(code_version=f'{base[:12]}+src.{digest[:16]}', source_sha256=digest,
                      executable_sha256=exe['sha256'], cpu_model=match.group(1).strip() if match else 'unknown',
                      host=exe['host'], slurm_job_id=self.job)
        submission = self.root/'submission.json'
        if submission.exists():
            saved = self.json(submission)
            check(saved['source_sha256'] == digest, 'Submission source identity mismatch')
            result['code_version'] = saved['code_version']
        for name in ['provenance.json', 'archive.json', 'slurm_accounting.txt']:
            path = self.root/name
            if path.exists():
                self.text(path)
        return result

    def stamp(self, row):
        row = {k: v for k, v in row.items() if k not in ('alpha_degrees', 'gamma')}
        for key, value in self.identity.items():
            check(key not in row or row[key] == value, f'Conflicting identity: {key}')
            row[key] = value
        return row

    def resources(self, prefix):
        self.json(self.data/f'{prefix}.command.json')
        value = self.json(self.data/f'{prefix}.resources.json')
        check(value['returncode'] == 0 and value['peak_rss_kib'] > 0, f'Failed process: {prefix}')
        check(value['host'] == self.identity['host'] and str(value['slurm_job_id']) == self.job,
              f'Process identity mismatch: {prefix}')
        return value

    def write(self, filename, records, lines, coverage):
        # Validation finishes before creating/updating any derived outputs.
        reports = sorted(p for p in self.root.glob('report_archived*') if p.is_dir())
        check(len(reports) <= 1, 'Multiple report directories; choose one canonical report first')
        report = reports[0] if reports else self.root/('report_archived_'+datetime.now(ZoneInfo('America/New_York')).strftime('%Y%m%d'))
        post = self.root/'postprocessing'
        post.mkdir(exist_ok=True)
        report.mkdir(exist_ok=True)
        table = post/filename
        with table.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
        # Maintain one current report and source manifest; no historical report copies.
        meta = report/'inputs.json'
        metadata = dict(schema_version=2, kind=self.kind, job_id=self.job, **coverage,
                        collected_utc=datetime.now(timezone.utc).isoformat(), identity=self.identity,
                        input_path_base='experiment_root',
                        inputs=[dict(path=str(p.relative_to(self.root)), sha256=sha(p)) for p in sorted(self.inputs)],
                        outputs=[dict(path=str(table.relative_to(self.root)), sha256=sha(table))],
                        processor_sha256=sha(Path(__file__)),
                        note='Raw inputs are referenced, never copied or modified. Historical missing metadata is not reconstructed beyond recorded evidence.')
        meta.write_text(json.dumps(metadata, indent=2, ensure_ascii=False)+'\n')
        links = [f'[后处理 CSV](../postprocessing/{filename})',
                 f'[原始数据](../job_{self.job}/data/)；[设备与构建环境](../job_{self.job}/environment.txt)。']
        environment = [f"CPU：{self.identity['cpu_model']}；主机：{self.identity['host']}；任务：{self.job}。", f"代码版本：`{self.identity['code_version']}`；二进制 SHA256：`{self.identity['executable_sha256']}`。"]
        text = '\n'.join(lines+['', *environment, '', *links, '', RESIDUAL_NOTE, '',
                 '数据目录与后处理文件不重复复制；输入路径相对于实验目录，完整输入哈希见 inputs.json。', ''])
        (report/'REPORT.md').write_text(text)
        (report/'SUCCESS').write_text(f'Validated requested {self.kind} coverage; full_formal_set={coverage["full_formal_set"]}\n')
        print(table)
        print(report/'REPORT.md')
        return report


def accuracy(root):
    run = Run(root, 'accuracy')
    records = run.csv(run.data/'convergence.csv')
    check(bool(records), 'Empty convergence table')
    records.sort(key=lambda r: (int(r['nx']), {'dg': 0, 'hdg': 1}[r['method']]))
    grids = sorted({int(r['nx']) for r in records})
    check([(int(r['nx']), r['method']) for r in records] == [(n, m) for n in grids for m in ('dg', 'hdg')], 'Missing/duplicate grid pairs')
    protocol = ('k', 'quad', 'cfl', 'integrator', 'linear', 'preconditioner', 'ordering', 'ntol', 'ltol',
                'safety', 'final_time', 'lx', 'ly')
    previous = {}
    raw = {}
    for n in grids:
        raw[n] = run.csv(run.data/f'n{n}_summary.csv')
        difference = run.csv(run.data/f'n{n}_difference.csv')
        check(len(difference) == 1, 'Missing difference record')
        check(float(difference[0]['max_nodal']) <= 1e-8, 'DG/HDG discrepancy')
        run.resources(f'n{n}')
    output = []
    for r in records:
        n = int(r['nx'])
        check(int(r['ny']) == n and int(r['fixed_newton_steps']) == 0, 'Expected converged square-grid run')
        check(all(r[key] == records[0][key] for key in protocol), 'Incompatible convergence protocol')
        matches = [x for x in raw[n] if x['method'] == r['method']]
        check(len(matches) == 1, 'Missing method summary')
        for component in COMPONENTS:
            key = r['method'], component
            error = float(r[f'{component}_L2'])
            check(math.isfinite(error) and error > 0, 'Invalid error')
            check(matches[0][f'{component}_L2'] == r[f'{component}_L2'], 'Raw/aggregate error mismatch')
            order = math.log(previous[key][1]/error)/math.log(n/previous[key][0]) if key in previous else None
            check(order is None or close(r[f'{component}_order'], order), 'Incorrect stored order')
            r[f'{component}_order'] = '' if order is None else order
            previous[key] = n, error
        first = {key: r[key] for key in ['nx', 'ny', 'method', 'k']}
        for component in COMPONENTS:
            for suffix in ['L2', 'order']:
                key = f'{component}_{suffix}'; first[key] = r[key]
        output.append(run.stamp(dict(first, **{k: v for k, v in r.items() if k not in first})))
    formal = grids == [32, 64, 128, 256] and int(records[0]['k']) == 2 and records[0]['integrator'] == 'ssp-dirk3' and close(records[0]['cfl'], 2) and close(records[0]['final_time'], .2)
    lines = ['# 2D 精度与收敛阶报告', '', f'作业 {run.job}；网格 {grids}；完整正式覆盖：{formal}。', '',
             '| N | 方法 | rho L2 | 阶 | mx L2 | 阶 | my L2 | 阶 | E L2 | 阶 |',
             '|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in output:
        values = [r['nx'], r['method']]
        for c in COMPONENTS:
            order = r[f'{c}_order']; values += [f"{float(r[f'{c}_L2']):.8e}", f'{order:.4f}' if order != '' else '—']
        lines.append('| '+' | '.join(values)+' |')
    lines += ['', f"DG/HDG 最大节点差：{max(float(r['difference_max_nodal']) for r in records):.4e}。"]
    return run.write('convergence_summary.csv', output, lines, dict(grids=grids, full_formal_set=formal))


def speed(root):
    run = Run(root, 'speed')
    timings = run.csv(run.data/'timings.csv')
    ratios = run.csv(run.data/'speedup.csv')
    check(bool(timings) and bool(ratios), 'Empty timing tables')
    degrees = sorted({int(r['k']) for r in timings}, reverse=True)
    check(sorted((int(r['k']), r['method']) for r in timings) == sorted((k, m) for k in degrees for m in ('dg', 'hdg')), 'Missing/duplicate timing pair')
    check(sorted(int(r['k']) for r in ratios) == sorted(degrees), 'Missing speedup pair')
    means, pair = {}, {}
    for r in timings:
        k, method = int(r['k']), r['method']
        prefix = f'k{k}_{method}'
        summary = run.csv(run.data/f'{prefix}_summary.csv')
        check(len(summary) == 1, 'Missing raw summary')
        check(all(r[key] == value for key, value in summary[0].items()), 'Raw/aggregate summary mismatch')
        check(int(r['steps']) == 1 and int(r['newton_iterations']) == 4 and int(r['fixed_newton_steps']) == 4 and r['integrator'] == 'be', 'Expected four BE updates')
        check(int(r['warmup_updates']) == 1 and int(r['measured_updates']) == 3, 'Incorrect warmup count')
        updates = run.csv(run.data/f'{prefix}_iterations.csv')
        check(len(updates) == 4 and [x['measured'] for x in updates] == ['0', '1', '1', '1'], 'Incorrect timing samples')
        mean = sum(float(x['iteration_s']) for x in updates[1:])/3
        check(mean > 0 and close(r['mean_newton_s'], mean), 'Incorrect Newton mean')
        metrics = run.resources(prefix)
        for key in ['wall_s', 'user_cpu_s', 'system_cpu_s', 'peak_rss_kib']:
            check(close(r[key], metrics[key]), f'Incorrect resource value: {key}')
        run.stamp(r)
        means[k, method] = mean
        pair[k, method] = r
    output = []
    for k in degrees:
        ratio = next(r for r in ratios if int(r['k']) == k)
        check(close(ratio['dg_mean_s'], means[k, 'dg']) and close(ratio['hdg_mean_s'], means[k, 'hdg']) and close(ratio['dg_over_hdg'], means[k, 'dg']/means[k, 'hdg']), 'Incorrect speedup')
        a, b = pair[k, 'dg'], pair[k, 'hdg']
        params = ['nx', 'ny', 'k', 'quad', 'cfl', 'integrator', 'linear', 'preconditioner', 'ordering', 'ntol', 'ltol', 'safety', 'fixed_newton_steps', 'warmup_updates']
        check(all(a[key] == b[key] for key in params), 'Incompatible paired timing protocol')
        row = dict(k=k, nx=a['nx'], ny=a['ny'], dg_mean_s=means[k, 'dg'], hdg_mean_s=means[k, 'hdg'], dg_over_hdg=means[k, 'dg']/means[k, 'hdg'])
        for method, r in [('dg', a), ('hdg', b)]:
            for key in ['wall_s', 'time_loop_s', 'user_cpu_s', 'system_cpu_s', 'krylov_iterations', 'max_linear_relative_residual', 'max_dg_linear_relative_residual']:
                row[f'{method}_{key}'] = r[key]
            row[f'{method}_peak_rss_gib'] = int(r['peak_rss_kib'])/1024**2
        row.update({key: a[key] for key in params if key not in row})
        output.append(run.stamp(row))
    formal = degrees == list(range(8, -1, -1)) and all(int(r['nx']) == 256 and int(r['ny']) == 256 and close(r['cfl'], 5) for r in timings)
    lines = ['# 2D 速度与内存报告', '', f'作业 {run.job}；阶数 {degrees}；完整正式覆盖：{formal}。', '',
             '均值仅统计后三次更新，Krylov总数和time_loop_s包含预热；wall_s是进程总时长。', '',
             '| k | DG s/更新 | HDG s/更新 | DG/HDG | DG RSS GiB | HDG RSS GiB |',
             '|---:|---:|---:|---:|---:|---:|']
    for r in output:
        lines.append(f"| {r['k']} | {r['dg_mean_s']:.6f} | {r['hdg_mean_s']:.6f} | {r['dg_over_hdg']:.4f} | {r['dg_peak_rss_gib']:.4f} | {r['hdg_peak_rss_gib']:.4f} |")
    lines += ['', '固定四次Newton更新不代表非线性收敛；timing-only误差字段nan不作为精度结论。']
    return run.write('speed_summary.csv', output, lines, dict(degrees=degrees, full_formal_set=formal))
