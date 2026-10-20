#!/usr/bin/env python3
"""512 cells: exactly four Newton updates, discard the first, average three."""
import argparse
import math
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
from experiment import prepare_experiment, rows, run, write_csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', required=True)
    parser.add_argument('--output', help='Output directory inside running/; default: unique local run')
    parser.add_argument('--degrees', nargs='+', type=int, default=list(range(9)))
    args = parser.parse_args()
    if len(set(args.degrees)) != len(args.degrees) or any(k not in range(9) for k in args.degrees):
        parser.error('degrees must be unique values in 0..8')
    output, exe = prepare_experiment(args.output, args.exe, 'speed', '1d')
    records, ratios = [], []
    for degree in args.degrees:
        pair = {}
        for method in ('dg', 'hdg'):
            prefix = output/f'k{degree}_{method}'
            resources = run([exe, '--cells', 512, '--k', degree, '--method', method,
                '--integrator', 'be', '--steps', 1, '--newton-steps', 4,
                '--warmup-updates', 1, '--cfl', 5, '--precond', 'block-jacobi',
                '--ordering', 'natural', '--ntol', '1e-12', '--ltol', '1e-12',
                '--restart', 20, '--max-linear', 1500, '--safety', 1.1,
                '--timing-only', '--prefix', prefix], prefix)
            summary = rows(f'{prefix}_summary.csv')[0]
            iterations = rows(f'{prefix}_iterations.csv')
            if (len(iterations) != 4 or int(summary['newton_iterations']) != 4
                    or int(summary['measured_updates']) != 3
                    or [r['measured'] for r in iterations] != ['0', '1', '1', '1']
                    or summary['integrator'] != 'be' or int(summary['steps']) != 1):
                raise RuntimeError('Expected four updates with exactly one warmup')
            mean = sum(float(r['iteration_s']) for r in iterations[1:])/3
            if not math.isfinite(mean) or mean <= 0 or not math.isclose(
                    mean, float(summary['mean_newton_s']), rel_tol=1e-13):
                raise RuntimeError('Incorrect measured Newton mean')
            summary.update(resources)
            summary['warmup_s'] = iterations[0]['iteration_s']
            records.append(summary)
            pair[method] = mean
            write_csv(output/'timings.csv', records)
        ratio = pair['dg']/pair['hdg']
        ratios.append(dict(k=degree, cells=512, cfl=5, integrator='be',
                           updates=4, warmup=1, measured_updates=3,
                           dg_mean_s=pair['dg'], hdg_mean_s=pair['hdg'],
                           dg_over_hdg=ratio, target_met=(ratio >= 8) if degree == 8 else ''))
        write_csv(output/'speedup.csv', ratios)
        print(f'k={degree}: DG/HDG={ratio:.4f}', flush=True)
    (output/'SUCCESS').write_text('All requested sequential DG/HDG timing pairs completed.\n')


if __name__ == '__main__':
    main()
