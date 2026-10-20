#!/usr/bin/env python3
"""Paired 2D Q2 SSP-DIRK3, CFL=2, T=0.2 convergence sweep."""
import argparse
import math
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from experiment import prepare_experiment,rows,run,write_csv

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe',required=True)
    parser.add_argument('--output', help='Output directory inside running/; default: unique local run')
    parser.add_argument('--grids',nargs='+',type=int,default=[32,64,128,256])
    args=parser.parse_args()
    if any(n<3 for n in args.grids) or args.grids!=sorted(set(args.grids)):
        parser.error('Grids must be distinct and increasing, N>=3')
    output,exe=prepare_experiment(args.output, args.exe, 'accuracy', '2d')
    records=[]; previous={}
    for n in args.grids:
        prefix=output/f'n{n}'
        print(f'Start N={n}, Q2, DIRK3, CFL=2, T=0.2',flush=True)
        run([exe,'--nx',n,'--ny',n,'--k',2,'--integrator','ssp-dirk3','--cfl',2,'--T',.2,
             '--method','both','--precond','block-jacobi',
             '--ordering','natural','--ntol','1e-12','--ltol','1e-12','--restart',20,
             '--max-linear',1500,'--max-newton',30,'--safety',1.1,'--alpha-deg',45,
             '--no-save-solution','--prefix',prefix],prefix)
        difference=rows(f'{prefix}_difference.csv')[0]
        summaries=rows(f'{prefix}_summary.csv')
        if [r['method'] for r in summaries]!=['dg','hdg']:raise RuntimeError('Missing paired results')
        for data in summaries:
            if float(data['final_time'])!=.2 or data['integrator']!='ssp-dirk3' or int(data['fixed_newton_steps'])!=0:
                raise RuntimeError('Incorrect convergence protocol')
            method=data['method']
            for component in ['rho','momentum_x','momentum_y','energy']:
                error=float(data[f'{component}_L2'])
                if not math.isfinite(error) or error<=0:raise RuntimeError('Invalid error')
                key=(method,component)
                data[f'{component}_order']=math.log(previous[key][1]/error)/math.log(n/previous[key][0]) if key in previous else ''
                previous[key]=(n,error)
                data[f'difference_{component}_L2']=difference[f'{component}_L2']
            data['difference_joint_L2']=difference['joint_L2']
            data['difference_max_nodal']=difference['max_nodal']
            records.append(data)
        write_csv(output/'convergence.csv',records)
        print(f'N={n} completed',flush=True)
    (output/'SUCCESS').write_text('All requested paired 2D grids completed.\n')
if __name__=='__main__':main()
