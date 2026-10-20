#!/usr/bin/env python3
"""Production DIRK3 and BE DG/HDG comparisons with exact block Jacobi."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from experiment import require_debug_build,prepare,rows,run

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe',required=True);parser.add_argument('--output',required=True)
    args=parser.parse_args();output,exe=prepare(args.output,require_debug_build(args.exe))
    for name,extra in [('dirk3',['--T',.02]),('be',['--integrator','be','--T',.02]),
                        ('entropy',['--testcase','wave','--T',.02])]:
        prefix=output/name
        run([exe,'--nx',4,'--ny',4,'--k',2,'--cfl',2,'--method','both','--check-update',
             '--ntol','1e-12','--ltol','1e-12','--no-save-solution',*extra,'--prefix',prefix],prefix)
        diff=rows(f'{prefix}_difference.csv')[0]
        if float(diff['max_nodal'])>1e-9:raise RuntimeError(f'{name}: DG/HDG discrepancy')
        for summary in rows(f'{prefix}_summary.csv'):
            if float(summary['max_conservation_drift'])>1e-10:raise RuntimeError(f'{name}: conservation')
        print(f'{name} passed',flush=True)
    (output/'SUCCESS').write_text('All production-path 2D checks passed.\n')
if __name__=='__main__':main()
