#!/usr/bin/env python3
"""2D 256x256 paired timings, k=8 down to 0; four updates, first warmup."""
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
    parser.add_argument('--degrees',type=int,nargs='+',default=list(range(8,-1,-1)))
    parser.add_argument('--cells',type=int,default=256)
    args=parser.parse_args()
    if args.cells<3 or len(set(args.degrees))!=len(args.degrees) or any(k not in range(9) for k in args.degrees):
        parser.error('Invalid mesh or degree list')
    output,exe=prepare_experiment(args.output, args.exe, 'speed', '2d')
    records=[];ratios=[]
    for k in args.degrees:
        pair={}
        for method in ['dg','hdg']:
            prefix=output/f'k{k}_{method}'
            print(f'Start {method}: N={args.cells}x{args.cells}, k={k}, exact block Jacobi',flush=True)
            resources=run([exe,'--nx',args.cells,'--ny',args.cells,'--k',k,'--method',method,
                '--integrator','be','--steps',1,'--newton-steps',4,'--warmup-updates',1,
                '--cfl',5,'--precond','block-jacobi',
                '--ordering','natural','--ntol','1e-12','--ltol','1e-12','--restart',20,
                '--max-linear',1500,'--safety',1.1,'--alpha-deg',45,'--timing-only','--prefix',prefix],prefix)
            summary=rows(f'{prefix}_summary.csv')[0]; updates=rows(f'{prefix}_iterations.csv')
            if len(updates)!=4 or [r['measured'] for r in updates]!=['0','1','1','1'] or int(summary['measured_updates'])!=3:
                raise RuntimeError('Expected exactly four updates, first warmup')
            mean=sum(float(r['iteration_s']) for r in updates[1:])/3
            if not math.isfinite(mean) or mean<=0 or not math.isclose(mean,float(summary['mean_newton_s']),rel_tol=1e-13):
                raise RuntimeError('Incorrect mean Newton time')
            summary.update(resources);summary['warmup_s']=updates[0]['iteration_s'];records.append(summary)
            pair[method]=mean;write_csv(output/'timings.csv',records)
            print(f'Completed {method} k={k}: mean={mean:.8g} s, RSS={resources["peak_rss_kib"]/1024/1024:.4f} GiB',flush=True)
        ratios.append(dict(k=k,nx=args.cells,ny=args.cells,preconditioner='block-jacobi',
            dg_jacobi_local='dense',hdg_jacobi_local='face-dense',integrator='be',cfl=5,
            updates=4,warmup=1,measured_updates=3,dg_mean_s=pair['dg'],hdg_mean_s=pair['hdg'],dg_over_hdg=pair['dg']/pair['hdg']))
        write_csv(output/'speedup.csv',ratios)
    (output/'SUCCESS').write_text('All requested sequential 2D timing pairs completed.\n')
if __name__=='__main__':main()
