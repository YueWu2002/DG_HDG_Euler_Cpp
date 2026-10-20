#!/usr/bin/env python3
"""Production DIRK smoke and converged one-step BE DG/HDG comparison."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
from experiment import prepare, rows, run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output, exe = prepare(args.output, args.exe)
    for name, options in (
        ('dirk_smoke', ['--cells', 16, '--k', 2, '--integrator', 'ssp-dirk3', '--T', .02]),
        ('be_update', ['--cells', 32, '--k', 3, '--integrator', 'be', '--steps', 1])):
        prefix = output/name
        run([exe, *options, '--method', 'both', '--cfl', 2, '--precond', 'block-jacobi',
             '--ntol', '1e-12', '--ltol', '1e-12', '--check-update', '--prefix', prefix], prefix)
        difference = rows(f'{prefix}_difference.csv')[0]
        if float(difference['max_nodal']) > 1e-9:
            raise RuntimeError(f'{name}: DG/HDG discrepancy')
        for data in rows(f'{prefix}_summary.csv'):
            if float(data['max_conservation_drift']) > 1e-10:
                raise RuntimeError(f'{name}: conservation failure')
        print(f'{name} passed')
    (output/'SUCCESS').write_text('Both production-path checks passed.\n')


if __name__ == '__main__':
    main()
