#!/usr/bin/env python3
"""Complete paired DIRK3/CFL=2 one-dimensional convergence sweep."""
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
    args = parser.parse_args()
    output, exe = prepare_experiment(args.output, args.exe, 'accuracy', '1d')
    records, previous = [], {}
    for cells in (32, 64, 128, 256, 512):
        prefix = output/f'n{cells}'
        run([exe, '--cells', cells, '--k', 2, '--integrator', 'ssp-dirk3',
             '--cfl', 2, '--T', .2, '--method', 'both', '--precond', 'block-jacobi',
             '--ordering', 'natural', '--ntol', '1e-12', '--ltol', '1e-12',
             '--restart', 20, '--max-linear', 1500, '--max-newton', 30,
             '--safety', 1.1, '--prefix', prefix], prefix)
        difference = rows(f'{prefix}_difference.csv')[0]
        summaries = rows(f'{prefix}_summary.csv')
        if [r['method'] for r in summaries] != ['dg', 'hdg']:
            raise RuntimeError('Missing paired solver results')
        for data in summaries:
            method = data['method']
            if (float(data['final_time']) != .2 or data['integrator'] != 'ssp-dirk3'
                    or int(data['fixed_newton_steps']) != 0):
                raise RuntimeError('Unexpected convergence protocol')
            for component in ('rho', 'momentum', 'energy'):
                error = float(data[f'{component}_L2'])
                if not math.isfinite(error) or error <= 0:
                    raise RuntimeError('Nonfinite/zero convergence error')
                key = method, component
                data[f'{component}_order'] = (math.log2(previous[key]/error)
                                              if key in previous else '')
                previous[key] = error
                data[f'difference_{component}_L2'] = difference[f'{component}_L2']
            data['difference_joint_L2'] = difference['joint_L2']
            data['difference_max_nodal'] = difference['max_nodal']
            records.append(data)
        write_csv(output/'convergence.csv', records)
        print(f'N={cells} completed', flush=True)
    (output/'SUCCESS').write_text('All five paired, fully converged DIRK3 grids completed.\n')


if __name__ == '__main__':
    main()
