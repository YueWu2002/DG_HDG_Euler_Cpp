#!/usr/bin/env python3
"""Freeze a new 2D speed sweep; build once on its selected compute node."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
from speed_campaign import rows, save, sha


def prepare(node, cells, factor, minimum, restart=20, build_time="02:00:00"):
    project=Path(__file__).resolve().parents[1]
    history=project/'speed_test/results/2d_submission_20261019T205542_KhFC6z/job_7635952/data/timings.csv'
    root=Path(tempfile.mkdtemp(prefix=f'2d_n{cells}_restart{restart}_{datetime.now():%Y%m%dT%H%M%S}_',dir=project/'speed_test/running'))
    for name in ['source','control','baseline','bin']:(root/name).mkdir()
    shutil.copy2(history,root/'baseline/timings.csv')
    history_rows=rows(history)
    reference=project/'MEMORY_ESTIMATE.md'
    reference_text=reference.read_text()
    section=reference_text.split('### 2D：256×256，GMRES(20)',1)[1].split('## 精度测试实测',1)[0]
    resource_rows={}
    for line in section.splitlines():
        parts=[x.strip() for x in line.split('|')]
        if len(parts)==7 and parts[1].isdigit():
            resource_rows[int(parts[1])]=[float(x) for x in parts[2:6]]
    if len(resource_rows)!=9:raise ValueError('Incomplete MEMORY_ESTIMATE reference table')
    (root/'resource_reference.md').write_text(section)
    scale=(cells/256)**2
    if factor<scale:raise ValueError('Resource factor must cover the cell-count scaling')
    cases=[]
    for k in range(8,-1,-1):
        for method in ['dg','hdg']:
            old=next(r for r in history_rows if int(r['k'])==k and r['method']==method)
            ref=resource_rows[k]
            historical_rss=ref[0 if method=='dg' else 1]
            historical_wall=ref[2 if method=='dg' else 3]
            dofs=cells*cells*4*(k+1)**2 if method=='dg' else 2*cells*cells*4*(k+1)
            extra_basis=max(0,restart-20)*dofs*8/1024**3
            estimated_memory=historical_rss*scale+extra_basis
            memory=max(minimum,math.ceil(estimated_memory*(factor/scale)))
            minutes=max(1,math.ceil(historical_wall*factor/60))
            duration=f'{minutes//60:02}:{minutes%60:02}:00'
            cases.append(dict(key=f'k{k}_{method}',k=k,method=method,historical_wall_s=historical_wall,
                              historical_peak_rss_kib=historical_rss*1024**2,resource_multiplier=factor,
                              extra_gmres_basis_gib=extra_basis,estimated_peak_memory_gib=estimated_memory,
                              estimated_wall_s=historical_wall*scale,
                              requested_memory_gib=memory,requested_time=duration))
    for name in ['speed_campaign.py','experiment.py','provenance.py']:
        shutil.copy2(project/'scripts'/name,root/'control'/name)
    (root/'control/job.sh').write_text('''#!/usr/bin/env bash
set -euo pipefail
module purge
module load intel-oneapi-compilers intel-oneapi-mkl eigen cmake
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
exec srun --export=ALL --ntasks=1 --cpus-per-task=1 --cpu-bind=cores python3 "$@"
''')
    (root/'control/postprocess_job.sh').write_text('''#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
exec srun --ntasks=1 --cpus-per-task=1 /usr/bin/python3 "$@"
''')
    files=['CMakeLists.txt','cmake','common','dg_hdg_1d','dg_hdg_2d','vendor','scripts','intel-build.sh','skills',
           'AGENTS.md','README.md','MATHEMATICS.md','DESIGN.md','VALIDATION.md','MEMORY_ESTIMATE.md','RESULTS.md',
           'code_test/2d','accuracy_test/2d','speed_test/2d']
    producer=subprocess.Popen(['tar','-C',str(project),'--exclude=__pycache__','--exclude=results','--exclude=running','-cf','-',*files],stdout=subprocess.PIPE)
    subprocess.run(['tar','-C',str(root/'source'),'-xf','-'],stdin=producer.stdout,check=True)
    producer.stdout.close()
    if producer.wait()!=0:raise RuntimeError('Snapshot creation failed')
    (root/'base_commit.txt').write_text(subprocess.check_output(['git','-C',str(project),'rev-parse','HEAD'],text=True))
    (root/'working_tree.patch').write_bytes(subprocess.check_output(['git','-C',str(project),'diff','--binary']))
    source=''.join(f'{sha(p)}  ./{p.relative_to(root/"source")}\n' for p in sorted((root/'source').rglob('*')) if p.is_file())
    (root/'source.sha256').write_text(source)
    digest=sha(root/'source.sha256');base=(root/'base_commit.txt').read_text().strip()
    driver=''.join(f'{sha(p)}  control/{p.name}\n' for p in sorted((root/'control').iterdir()))
    (root/'driver.sha256').write_text(driver)
    save(root/'submission.json',dict(schema_version=1,experiment='speed',dimension='2d',cells=cells,restart=restart,
         fresh_build=True,build_time=build_time,build_resource_note='16GiB; two configurations plus Debug checks; 2h initial budget with uncertain Debug runtime',expected_node=node,expected_cpu=None,executable_sha256=None,code_version=f'{base[:12]}+src.{digest[:16]}',
         source_sha256=digest,base_commit=base,driver_manifest_sha256=sha(root/'driver.sha256'),cases=cases,
         baseline_job_id='7635952',baseline_restart=20,baseline_timings_sha256=sha(root/'baseline/timings.csv'),
         minimum_memory_gib=minimum,resource_rule='MEMORY_ESTIMATE table; RSS times cells ratio plus extra restart basis, then margin; wall times factor, rounded up',
         resource_reference='MEMORY_ESTIMATE.md / 2D 256x256 GMRES20',resource_reference_sha256=sha(reference),
         time_uncertainty='Restart change may change iteration count and orthogonalization cost; no measured 512x512 GMRES30 baseline',
         retain_source_after_completion=False,source_directory='source',source_manifest='source.sha256',
         created_utc=datetime.now(timezone.utc).isoformat(),scheduling_note='Each (k,method) independently scheduled on one node after shared build; concurrency varies.'))
    print(root)
    for c in cases:print(c['key'],c['requested_memory_gib'],'GiB',c['requested_time'])
    return root


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--node',required=True)
    parser.add_argument('--cells',type=int,default=512)
    parser.add_argument('--restart',type=int,default=20)
    parser.add_argument('--build-time',default='02:00:00')
    parser.add_argument('--factor',type=float,default=4.4)
    parser.add_argument('--minimum-memory-gib',type=int,default=16)
    args=parser.parse_args()
    if args.cells<3 or args.factor<=0 or args.minimum_memory_gib<1 or args.restart<1:parser.error('Invalid resource or mesh request')
    prepare(args.node,args.cells,args.factor,args.minimum_memory_gib,args.restart,args.build_time)
