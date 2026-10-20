#!/usr/bin/env python3
"""Validate the full 2D experiment protocol and archive a standalone report."""
import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import shutil


def rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def check(condition, message):
    if not condition:
        raise ValueError(message)


def close(value, target):
    return math.isfinite(float(value)) and math.isclose(float(value), target, rel_tol=1e-12, abs_tol=1e-15)


def protocol(row, n, k, integrator, cfl, fixed_updates):
    check(int(row['nx']) == n and int(row['ny']) == n and int(row['k']) == k, 'Incorrect grid/degree')
    check(int(row['quad']) == k+1, 'Expected Gauss collocation')
    check(row['linear'] == 'gmres' and row['preconditioner'] == 'block-jacobi', 'Incorrect linear solver')
    check(row['ordering'] == 'natural' and row['integrator'] == integrator, 'Incorrect integration/ordering')
    check(close(row['cfl'], cfl) and close(row['ntol'], 1e-12) and close(row['ltol'], 1e-12), 'Incorrect tolerances/CFL')
    check(close(row['safety'], 1.1) and close(row['alpha_degrees'], 45) and close(row['gamma'], 3), 'Incorrect example')
    check(int(row['fixed_newton_steps']) == fixed_updates, 'Incorrect Newton protocol')
    expected_local = 'dense' if row['method'] == 'dg' else 'face-dense'
    check(row['jacobi_local'] == expected_local, 'Formal report requires exact block Jacobi')
    for name in ['max_linear_relative_residual', 'max_dg_linear_relative_residual', 'max_conservation_drift']:
        check(math.isfinite(float(row[name])) and float(row[name]) >= 0, f'Invalid {name}')


def resources(data, name):
    metrics = json.loads((data/f'{name}.resources.json').read_text())
    check(metrics['returncode'] == 0 and metrics['peak_rss_kib'] > 0, f'Failed resource record: {name}')
    return metrics


def collect(accuracy_job, speed_job, output):
    accuracy = accuracy_job/'data'
    speed = speed_job/'data'
    for data in [accuracy, speed]:
        check((data/'SUCCESS').is_file(), f'Experiment is incomplete: {data}')
    convergence = rows(accuracy/'convergence.csv')
    check([(int(r['nx']), r['method']) for r in convergence] ==
          [(n, method) for n in [32, 64, 128, 256] for method in ['dg', 'hdg']], 'Missing convergence pairs')
    previous = {}
    accuracy_resources = []
    for row in convergence:
        n = int(row['nx'])
        protocol(row, n, 2, 'ssp-dirk3', 2, 0)
        check(close(row['final_time'], .2), 'Incorrect final accuracy time')
        check(float(row['difference_max_nodal']) <= 1e-8, 'DG/HDG mismatch')
        for component in ['rho', 'momentum_x', 'momentum_y', 'energy']:
            error = float(row[f'{component}_L2'])
            check(math.isfinite(error) and error > 0, 'Invalid convergence error')
            key = (row['method'], component)
            if key in previous:
                expected = math.log2(previous[key]/error)
                check(close(row[f'{component}_order'], expected), 'Incorrect reported order')
            previous[key] = error
        if row['method'] == 'dg':
            accuracy_resources.append(resources(accuracy, f'n{n}'))

    timings = rows(speed/'timings.csv')
    ratios = rows(speed/'speedup.csv')
    check([(int(r['k']), r['method']) for r in timings] ==
          [(k, method) for k in range(8, -1, -1) for method in ['dg', 'hdg']], 'Missing timing pairs')
    check([int(r['k']) for r in ratios] == list(range(8, -1, -1)), 'Missing speedup rows')
    means = {}
    for row in timings:
        k = int(row['k'])
        protocol(row, 256, k, 'be', 5, 4)
        check(int(row['steps']) == 1 and int(row['newton_iterations']) == 4, 'Incorrect update count')
        check(int(row['warmup_updates']) == 1 and int(row['measured_updates']) == 3, 'Incorrect warmup protocol')
        name = f'k{k}_{row["method"]}'
        updates = rows(speed/f'{name}_iterations.csv')
        check(len(updates) == 4 and [r['measured'] for r in updates] == ['0', '1', '1', '1'], 'Incorrect timing samples')
        mean = sum(float(r['iteration_s']) for r in updates[1:])/3
        check(mean > 0 and close(row['mean_newton_s'], mean), 'Incorrect measured mean')
        means[(k, row['method'])] = mean
        metrics = resources(speed, name)
        check(row['host'] == metrics['host'] and row['slurm_job_id'] == metrics['slurm_job_id'], 'Incorrect process metadata')
    check(len({r['host'] for r in timings}) == 1 and len({r['slurm_job_id'] for r in timings}) == 1,
          'Timing pairs did not run on the same node/job')
    for row in ratios:
        k = int(row['k'])
        check(close(row['dg_mean_s'], means[(k, 'dg')]) and close(row['hdg_mean_s'], means[(k, 'hdg')]), 'Inconsistent timing tables')
        check(close(row['dg_over_hdg'], means[(k, 'dg')]/means[(k, 'hdg')]), 'Incorrect speedup')

    output.mkdir(parents=True, exist_ok=False)
    shutil.copytree(accuracy, output/'accuracy')
    shutil.copytree(speed, output/'speed')
    for label, root in [('accuracy', accuracy_job), ('speed', speed_job)]:
        provenance = output/f'{label}_provenance'
        provenance.mkdir()
        for name in ['environment.txt', 'ctest.log', 'code_test.log']:
            shutil.copy2(root/name, provenance/name)
        for name in ['source.sha256', 'base_commit.txt', 'job_id.txt']:
            shutil.copy2(root.parent/name, provenance/name)
    (output/'inputs.json').write_text(json.dumps(dict(
        collected_utc=datetime.now(timezone.utc).isoformat(),
        accuracy_job=str(accuracy_job), speed_job=str(speed_job)), indent=2)+'\n')
    lines = ['# 2D 完整实验报告', '',
        '已验收四组 DG/HDG 收敛网格及九组次数的完整速度样本。', '',
        '精度：Q2、DIRK3、CFL=2、T=0.2。速度：256²、BE、CFL=5，精确右块 Jacobi，四次 Newton 更新首个预热。', '',
        '## 收敛', '', '| N | 方法 | rho L2 | 阶 | mx L2 | 阶 | my L2 | 阶 | E L2 | 阶 |',
        '|---|---|---|---|---|---|---|---|---|---|']
    for row in convergence:
        values = [row['nx'], row['method']]
        for component in ['rho', 'momentum_x', 'momentum_y', 'energy']:
            order = row[f'{component}_order']
            values += [f'{float(row[f"{component}_L2"]):.8e}', f'{float(order):.4f}' if order else '—']
        lines.append('| '+' | '.join(values)+' |')
    lines += ['', f'DG/HDG 最大节点差：{max(float(r["difference_max_nodal"]) for r in convergence):.6e}。', '',
        '## 速度与内存', '', '| k | DG 均值/s | HDG 均值/s | DG/HDG | DG 峰值/GiB | HDG 峰值/GiB | DG/HDG Krylov总数 |',
        '|---|---|---|---|---|---|---|']
    for row in ratios:
        k = int(row['k'])
        pair = [r for r in timings if int(r['k']) == k]
        values = [str(k), f'{means[(k,"dg")]:.8g}', f'{means[(k,"hdg")]:.8g}', f'{float(row["dg_over_hdg"]):.4f}']
        values += [f'{int(r["peak_rss_kib"])/1024**2:.4f}' for r in pair]
        values += ['/'.join(r['krylov_iterations'] for r in pair)]
        lines.append('| '+' | '.join(values)+' |')
    lines += ['', 'Krylov 总数包含预热更新；均值仅包含后三次。', '',
        '[完整收敛 CSV](accuracy/convergence.csv)、[速度比 CSV](speed/speedup.csv)、[逐方法资源与时间](speed/timings.csv)。', '',
        '原始命令、每次更新、逐进程 RSS、源码与二进制哈希随本报告归档。',
        '上述数据使用精确单元/面块 LU，只说明本次机器与参数下的结果。', '']
    (output/'REPORT.md').write_text('\n'.join(lines))
    (output/'SUCCESS').write_text('Full 2D accuracy and speed protocols validated.\n')
    print(output/'REPORT.md')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--accuracy-job-root', type=Path, required=True)
    parser.add_argument('--speed-job-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        collect(args.accuracy_job_root.resolve(), args.speed_job_root.resolve(), args.output.resolve())
    except (ValueError, OSError) as error:
        parser.exit(1, f'Cannot collect complete 2D results: {error}\n')


if __name__ == '__main__':
    main()
