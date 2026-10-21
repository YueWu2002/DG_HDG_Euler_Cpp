#!/usr/bin/env python3
"""Independent single-method Slurm timings with dependent collection/archival."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import getpass
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

CPU = 'Intel(R) Xeon(R) Platinum 8268 CPU @ 2.90GHz'
EXCLUDED = {'alpha_degrees', 'gamma'}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def save(path, data):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False)+'\n')
    temporary.replace(path)


def rows(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def table(path, data):
    if not data:
        return
    with Path(path).open('w', newline='') as stream:
        fields = list(dict.fromkeys(k for row in data for k in row))
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(data)


def ensure(ok, message):
    if not ok:
        raise ValueError(message)


def equal(a, b):
    return math.isfinite(float(a)) and math.isclose(float(a), float(b), rel_tol=1e-11, abs_tol=1e-14)


def read_manifest(root):
    manifest = json.loads((root/'submission.json').read_text())
    if manifest.get('fresh_build') and (root/'build_identity.json').exists():
        manifest.update(json.loads((root/'build_identity.json').read_text()))
    return manifest


def baseline(dim, project):
    # Read-only historical paths: do not rename archived evidence. New copies use bin/solver.exe.
    if dim == '1d':
        root = project/'speed_test/results/1d_slurm_7622721'
        return dict(node='node1837', binary=project/'build-intel/dg_hdg_1d/dg_hdg',
                    expected='da18d8c6e359e544259dd911f8e43ae08f757ed6b5dc5d83a01c5ce54357f8f2',
                    timings=root/'timings_slurm_7622721.csv', baseline_job='7622721', source_manifest=None,
                    environment=None, build=project/'build-intel')
    root = project/'speed_test/results/2d_submission_20261019T205542_KhFC6z'
    return dict(node='node2301', binary=root/'job_7635952/build/dg_hdg_2d/dg_hdg_2d',
                expected='812e1075efc2854e89cfc02de8426e23ed757eb1abd497fc104ce2771e98b533',
                timings=root/'job_7635952/data/timings.csv', baseline_job='7635952',
                source_manifest=root/'source.sha256', environment=root/'job_7635952/environment.txt',
                build=root/'job_7635952/build')


def prepare(dim, restart, minimum, cells=512, node=None):
    project = Path(__file__).resolve().parents[1]
    b = baseline(dim, project)
    ensure(sha(b['binary']) == b['expected'], f'Baseline binary identity mismatch: {dim}')
    running = project/'speed_test/running'
    running.mkdir(exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix=f'{dim}_restart{restart}_{datetime.now():%Y%m%dT%H%M%S}_', dir=running))
    for name in ['bin', 'control', 'baseline']:
        (root/name).mkdir()
    shutil.copy2(b['binary'], root/'bin/solver.exe')
    shutil.copy2(b['timings'], root/'baseline/timings.csv')
    for name in ['speed_campaign.py', 'experiment.py', 'provenance.py']:
        shutil.copy2(project/'scripts'/name, root/'control'/name)
    for name in ['CMakeCache.txt', 'compile_commands.json']:
        if (b['build']/name).exists():
            shutil.copy2(b['build']/name, root/'baseline'/name)
    if b['environment']:
        shutil.copy2(b['environment'], root/'baseline/environment.txt')
    source_sha = None
    base_commit = None
    if b['source_manifest']:
        shutil.copy2(b['source_manifest'], root/'source.sha256')
        source_sha = sha(root/'source.sha256')
        commit_file = b['source_manifest'].parent/'base_commit.txt'
        base_commit = commit_file.read_text().strip()
        shutil.copy2(commit_file, root/'base_commit.txt')
    driver_manifest = ''.join(f'{sha(p)}  control/{p.name}\n' for p in sorted((root/'control').iterdir()))
    (root/'driver.sha256').write_text(driver_manifest)
    degrees = list(range(9)) if dim == '1d' else list(range(8, -1, -1))
    history = {(int(r['k']), r['method']): r for r in rows(b['timings'])}
    cases = []
    for k in degrees:
        for method in ['dg', 'hdg']:
            r = history[k, method]
            ensure(r['host'].split('.')[0] == b['node'], 'Historical node mismatch')
            n = (cells*3*(k+1) if method == 'dg' else cells*3) if dim == '1d' else (
                256**2*4*(k+1)**2 if method == 'dg' else 2*256**2*4*(k+1))
            offset = max(0, restart-20)*n*8/1024**3
            old_rss = int(r['peak_rss_kib'])/1024**2 * (max(1, cells/512) if dim == '1d' else 1)
            memory = max(minimum, math.ceil((old_rss+offset)*1.2+1))
            cases.append(dict(key=f'k{k}_{method}', k=k, method=method, historical_rss_gib=old_rss,
                              extra_basis_gib=offset, requested_memory_gib=memory))
    manifest = dict(schema_version=1, experiment='speed', dimension=dim, restart=restart, baseline_restart=20,
                    baseline_job_id=b['baseline_job'], baseline_timings_sha256=sha(root/'baseline/timings.csv'),
                    cells=cells, expected_node=node or b['node'], expected_cpu=CPU, executable_sha256=b['expected'],
                    code_version=f'{base_commit[:12]}+src.{source_sha[:16]}' if source_sha else f'binary.{b["expected"][:16]}',
                    source_sha256=source_sha, base_commit=base_commit,
                    source_identity_note='Reuses byte-identical baseline solver. 1D historical source identity is unavailable; binary identity is verified.',
                    driver_manifest_sha256=sha(root/'driver.sha256'), cases=cases,
                    created_utc=datetime.now(timezone.utc).isoformat(), minimum_memory_gib=minimum,
                    parameter_change=f'cells={cells}, restart={restart}; independent method tasks on one selected node',
                    scheduling_note='Concurrency and node co-tenancy differ from the sequential baseline; record contention in timing interpretation.',
                    retain_source_after_completion=False)
    save(root/'submission.json', manifest)
    lines = ['#!/usr/bin/env bash', 'set -euo pipefail',
             'module purge', 'module load intel-oneapi-compilers intel-oneapi-mkl eigen cmake',
             'export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1',
             'export VECLIB_MAXIMUM_THREADS=1',
             'exec srun --export=ALL --ntasks=1 --cpus-per-task=1 --cpu-bind=cores python3 "$@"']
    (root/'control/job.sh').write_text('\n'.join(lines)+'\n')
    (root/'control/postprocess_job.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec srun --ntasks=1 --cpus-per-task=1 /usr/bin/python3 \"$@\"\n')
    driver_manifest = ''.join(f'{sha(p)}  control/{p.name}\n' for p in sorted((root/'control').iterdir()) if p.is_file())
    (root/'driver.sha256').write_text(driver_manifest)
    manifest['driver_manifest_sha256'] = sha(root/'driver.sha256')
    save(root/'submission.json', manifest)
    return root


def sbatch(root, name, memory, mode, extra=(), dependency=None, node=None, duration='1-00:00:00', hold=False):
    logs = root/'scheduler'; logs.mkdir(exist_ok=True)
    command = ['sbatch', '--parsable', '--partition=batch', '--nodes=1', '--ntasks=1', '--cpus-per-task=1',
               f'--mem={memory}G', f'--time={duration}', '--job-name='+name, '--chdir='+str(root.parent.parent if mode=='archive' else root),
               '--output='+str(logs/'%x-%j.out'), '--error='+str(logs/'%x-%j.err'), '--no-requeue']
    if hold:
        command += ['--hold']
    if dependency:
        command += ['--dependency='+dependency, '--kill-on-invalid-dep=yes']
    if node:
        command += ['--nodelist='+node, '--constraint=intel']
    shell = root/'control'/('job.sh' if mode in ('case','build') else 'postprocess_job.sh')
    runner = root/'control'/('speed_campaign.py' if mode in ('case','build') else 'postprocess_campaign.py')
    if not runner.exists():
        runner = root/'control/speed_campaign.py'
    command += [str(shell), str(runner), mode, '--run-root', str(root), *extra]
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    job = result.stdout.strip().split(';')[0]
    ensure(job.isdigit(), f'Invalid sbatch response: {result.stdout}')
    return job


def submit(root):
    manifest = read_manifest(root)
    index = root/'submitted_jobs.csv'
    if index.exists():
        raise ValueError('Submission already started; do not duplicate jobs')
    submitted = []
    build_dependency = None
    if manifest.get('node_probe_job_id'):
        submitted.append(dict(case='node_probe',job_id=manifest['node_probe_job_id'],memory_gib=16,node=manifest['expected_node'],role='probe'))
        table(index,submitted)
    if manifest.get('fresh_build'):
        build_job = sbatch(root, f'eigen-2d-n{manifest.get("cells",256)}-r{manifest["restart"]}-build', 16, 'build', node=manifest['expected_node'], duration=manifest.get('build_time','00:30:00'))
        submitted.append(dict(case='build',job_id=build_job,memory_gib=16,node=manifest['expected_node'],role='build'))
        table(index, submitted)
        build_dependency = 'afterok:'+build_job
    # Persist each successful submission immediately, including partial failures.
    for case in manifest['cases']:
        job = sbatch(root, f'eigen-{manifest["dimension"]}-r{manifest["restart"]}-{case["key"]}',
                     case['requested_memory_gib'], 'case', ['--case', case['key']], node=manifest['expected_node'], dependency=build_dependency, duration=case.get('requested_time','1-00:00:00'))
        submitted.append(dict(case=case['key'], job_id=job, memory_gib=case['requested_memory_gib'], node=manifest['expected_node'], role='compute'))
        table(index, submitted)
    ids = ':'.join(r['job_id'] for r in submitted)
    collector = sbatch(root, f'eigen-{manifest["dimension"]}-r{manifest["restart"]}-collect',
                       manifest['minimum_memory_gib'], 'collect', dependency='afterany:'+ids, duration='00:10:00')
    submitted.append(dict(case='collect',job_id=collector,memory_gib=manifest['minimum_memory_gib'],node='',role='collect'))
    table(index, submitted)
    # Archive only after the collector exits; never move a compute/collector's live working directory.
    archive = sbatch(root, f'eigen-{manifest["dimension"]}-r{manifest["restart"]}-archive',
                     manifest['minimum_memory_gib'], 'archive', dependency='afterany:'+collector, duration='00:10:00')
    submitted.append(dict(case='archive',job_id=archive,memory_gib=manifest['minimum_memory_gib'],node='',role='archive'))
    table(index, submitted)
    save(root/'job_ids.json', dict(compute=[r['job_id'] for r in submitted if r['role']=='compute'], collect=collector, archive=archive))
    print(json.dumps(dict(run_root=str(root),compute_tasks=len(manifest['cases']),collect_job=collector,archive_job=archive)))



def build(root):
    from provenance import runtime
    folder=root/f"job_{os.environ['SLURM_JOB_ID']}";folder.mkdir()
    save(folder/'provenance.json',runtime())
    (folder/'environment.txt').write_text(subprocess.check_output(['lscpu'],text=True))
    (folder/'slurm_allocation.txt').write_text(subprocess.check_output(['scontrol','show','job',os.environ['SLURM_JOB_ID']],text=True))
    env=dict(os.environ,BUILD_DIR=str(folder/'build'),BUILD_JOBS='1',BUILD_TESTING='OFF',BUILD_TYPE='Release',BUILD_1D='OFF',BUILD_2D='ON')
    status=1
    try:
        info = runtime()
        manifest = read_manifest(root)
        ensure(info['host'].split('.')[0] == manifest['expected_node'], 'Build node mismatch')
        if manifest.get('required_cpu'):
            ensure(info['cpu_model'] == manifest['required_cpu'], 'Build CPU model mismatch')
        subprocess.run(['bash',str(root/'source/scripts/build_debug_checks.sh'),'2d',str(folder)],env=env,check=True)
        with (folder/'build.log').open('w') as stream:
            subprocess.run(['bash',str(root/'source/intel-build.sh')],env=env,stdout=stream,stderr=subprocess.STDOUT,check=True)
        exe=folder/'build/bin/dg_hdg_2d.exe'
        ensure(runtime(exe)['cpu_model'] == read_manifest(root).get('required_cpu', runtime(exe)['cpu_model']), 'Build CPU model mismatch')
        shutil.copy2(exe,root/'bin/solver.exe')
        info=runtime(exe)
        ensure(info['host'].split('.')[0]==read_manifest(root)['expected_node'],'Build node mismatch')
        save(folder/'binary_provenance.json',info)
        with (folder/'environment.txt').open('a') as stream:
            subprocess.run(['bash','-c','icpx --version; cmake --version; module list'],stdout=stream,stderr=subprocess.STDOUT)
        save(root/'build_identity.json',dict(executable_sha256=sha(exe),expected_cpu=info['cpu_model'],build_job_id=os.environ['SLURM_JOB_ID'],build_directory=str((folder/'build').relative_to(root))))
        status=0
    finally:
        (folder/'exit_code.txt').write_text(str(status)+'\n')


def case(root, key):
    from experiment import prepare as prepare_data, run as run_command, write_csv
    from provenance import runtime
    manifest = read_manifest(root)
    spec = next(c for c in manifest['cases'] if c['key'] == key)
    job = os.environ['SLURM_JOB_ID']
    folder = root/f'job_{job}'; folder.mkdir()
    save(folder/'case.json', spec)
    try:
        env = runtime(root/'bin/solver.exe')
        save(folder/'provenance.json', env)
        (folder/'environment.txt').write_text(subprocess.check_output(['lscpu'], text=True))
        (folder/'slurm_allocation.txt').write_text(subprocess.check_output(['scontrol','show','job',job], text=True))
        ensure(env['cpu_model'] == manifest['expected_cpu'], 'Unexpected CPU model')
        ensure(env['host'].split('.')[0] == manifest['expected_node'], 'Unexpected execution node')
        ensure(sha(root/'bin/solver.exe') == manifest['executable_sha256'], 'Binary hash mismatch')
        output, exe = prepare_data(folder/'data', root/'bin/solver.exe')
        prefix = output/key
        command = [exe, '--cells',manifest.get('cells',512)] if manifest['dimension']=='1d' else [exe,'--nx',manifest.get('cells',256),'--ny',manifest.get('cells',256)]
        command += ['--k',spec['k'],'--method',spec['method'],'--integrator','be','--steps',1,
                    '--newton-steps',4,'--warmup-updates',1,'--cfl',5,'--precond','block-jacobi',
                    '--ordering','natural','--ntol','1e-12','--ltol','1e-12','--restart',manifest['restart'],
                    '--max-linear',1500,'--safety',1.1,'--timing-only','--prefix',prefix]
        resources = run_command(command, prefix)
        summary = rows(str(prefix)+'_summary.csv')[0]
        updates = rows(str(prefix)+'_iterations.csv')
        validate(summary, updates)
        record = {k:v for k,v in summary.items() if k not in EXCLUDED}
        record.update(resources);record['restart']=manifest['restart']
        write_csv(output/'timings.csv', [record])
        (output/'SUCCESS').write_text('One method, four updates; first warmup, three measured.\n')
        (folder/'exit_code.txt').write_text('0\n')
        print(key,summary['mean_newton_s'],flush=True)
    except BaseException:
        (folder/'exit_code.txt').write_text('1\n')
        raise


def validate(summary, updates):
    ensure(len(updates)==4 and [r['measured'] for r in updates]==['0','1','1','1'], 'Incorrect update samples')
    ensure(int(summary['measured_updates'])==3 and int(summary['newton_iterations'])==4 and int(summary['steps'])==1, 'Incorrect summary counts')
    mean = sum(float(r['iteration_s']) for r in updates[1:])/3
    ensure(mean>0 and equal(mean,summary['mean_newton_s']), 'Incorrect mean time')
    ensure(summary['integrator']=='be' and summary['preconditioner']=='block-jacobi' and equal(summary['cfl'],5), 'Incorrect protocol')
    return mean


def accounting(root):
    jobs = [r['job_id'] for r in rows(root/'submitted_jobs.csv') if r['role']!='archive']
    cmd=['sacct','-n','-P','--jobs',','.join(jobs),'--format=JobIDRaw,JobName,State,ExitCode,Elapsed,NodeList,MaxRSS']
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    (root/'slurm_accounting.txt').write_text(result.stdout)
    return result.stdout


def collect(root):
    manifest = read_manifest(root)
    submitted = rows(root/'submitted_jobs.csv')
    baseline_rows = {(int(r['k']),r['method']):r for r in rows(root/'baseline/timings.csv')}
    inputs = {root/'submission.json',root/'baseline/timings.csv',root/'submitted_jobs.csv'}
    records=[]; failures=[]
    for entry in submitted:
        if entry['role']!='compute':continue
        job=entry['job_id'];folder=root/f'job_{job}';data=folder/'data';key=entry['case']
        try:
            ensure((folder/'exit_code.txt').read_text().strip()=='0' and (data/'SUCCESS').exists(), 'Case failed/incomplete')
            summary=rows(data/f'{key}_summary.csv')[0]
            updates=rows(data/f'{key}_iterations.csv')
            validate(summary,updates)
            resources=json.loads((data/f'{key}.resources.json').read_text())
            provenance=json.loads((folder/'provenance.json').read_text())
            ensure(resources['returncode']==0 and provenance['executable_sha256']==manifest['executable_sha256'], 'Process/binary mismatch')
            ensure(provenance['cpu_model']==manifest['expected_cpu'] and resources['host'].split('.')[0]==manifest['expected_node'], 'CPU/node mismatch')
            command=json.loads((data/f'{key}.command.json').read_text())
            ensure(command[command.index('--restart')+1]==str(manifest['restart']), 'Incorrect restart')
            ensure('--alpha-deg' not in command, 'Unexpected example CLI parameter')
            old=baseline_rows[int(summary['k']),summary['method']]
            for param in ['k','quad','cfl','integrator','linear','preconditioner','ordering','ntol','ltol','safety','fixed_newton_steps','warmup_updates','steps']:
                ensure(summary[param]==old[param], f'Protocol changed: {param}')
            for param in (['cells'] if manifest['dimension']=='1d' else ['nx','ny']):
                ensure(int(summary[param]) == (manifest.get('cells',512) if param=='cells' else (manifest.get('cells',256) if manifest.get('fresh_build') else int(old[param]))), f'Incorrect mesh: {param}')
            record={k:v for k,v in summary.items() if k not in EXCLUDED}
            record.update(resources)
            record.update(restart=manifest['restart'],code_version=manifest['code_version'],source_sha256=manifest['source_sha256'],
                          executable_sha256=manifest['executable_sha256'],cpu_model=manifest['expected_cpu'],baseline_restart=20,
                          baseline_mean_newton_s=old['mean_newton_s'] if not manifest.get('fresh_build') and summary.get('cells')==old.get('cells') else '',old_over_new=float(old['mean_newton_s'])/float(summary['mean_newton_s']) if not manifest.get('fresh_build') and summary.get('cells')==old.get('cells') else '')
            records.append(record)
            inputs.update(p for p in data.iterdir() if p.suffix in ['.csv','.json'])
            inputs.add(folder/'provenance.json')
        except Exception as error:
            failures.append(dict(case=key,job_id=job,error=str(error)))
    post=root/'postprocessing';post.mkdir(exist_ok=True)
    report=root/'report_archived';report.mkdir(exist_ok=True)
    # Retain valid partial results if a task fails; never fill missing values with zero.
    table(post/'timings.csv',records)
    ratios=[]
    for k in sorted({int(r['k']) for r in records},reverse=True):
        pair={r['method']:r for r in records if int(r['k'])==k}
        if set(pair)!= {'dg','hdg'}:continue
        a,b=pair['dg'],pair['hdg']
        row=dict(k=k,restart=manifest['restart'],dg_mean_s=a['mean_newton_s'],hdg_mean_s=b['mean_newton_s'],
                 dg_over_hdg=float(a['mean_newton_s'])/float(b['mean_newton_s']),
                 code_version=manifest['code_version'],source_sha256=manifest['source_sha256'],
                 executable_sha256=manifest['executable_sha256'],cpu_model=manifest['expected_cpu'],node=manifest['expected_node'])
        for method,x in pair.items():
            for field in ['wall_s','time_loop_s','krylov_iterations','max_linear_relative_residual','max_dg_linear_relative_residual','slurm_job_id','old_over_new']:
                row[f'{method}_{field}']=x[field]
            row[f'{method}_peak_rss_gib']=int(x['peak_rss_kib'])/1024**2
        ratios.append(row)
    table(post/'speed_summary.csv',ratios)
    save(post/'failures.json',failures)
    complete=not failures and len(records)==18 and len(ratios)==9
    lines=[f'# {manifest["dimension"]} restart={manifest["restart"]} 速度报告','',
           f'完整覆盖：{complete}；原节点：{manifest["expected_node"]}；CPU：{manifest["expected_cpu"]}。',
           ('由本批不可变源码在指定计算节点编译；每方法四次Newton更新，首个预热、后三次平均。' if manifest.get('fresh_build') else '复用restart=20基线的同一二进制；每方法四次Newton更新，首个预热、后三次平均。'),
           '按用户要求分任务并发；同节点的资源争用与历史顺序执行不同，不能把全部耗时变化只归因于restart。','',
           '| k | DG s/更新 | HDG s/更新 | DG/HDG | DG旧/新 | HDG旧/新 |','|---:|---:|---:|---:|---:|---:|']
    for row in ratios:
        lines.append(f"| {row['k']} | {float(row['dg_mean_s']):.6g} | {float(row['hdg_mean_s']):.6g} | {row['dg_over_hdg']:.4f} | {row['dg_old_over_new']} | {row['hdg_old_over_new']} |")
    lines+=['','[后处理CSV](../postprocessing/speed_summary.csv)；[逐方法时间与资源](../postprocessing/timings.csv)；[失败记录](../postprocessing/failures.json)。',
            '原始输出保存在各job_<id>/data/；Slurm输出与错误保存在scheduler/，不复制原始数据。',
            'HDG沿用原Schur Krylov残差与恢复后DG验收规则，未改变容差或正交化。','']
    (report/'REPORT.md').write_text('\n'.join(lines))
    save(report/'inputs.json',dict(schema_version=2,dimension=manifest['dimension'],restart=manifest['restart'],full_formal_set=complete,
         inputs=[dict(path=str(p.relative_to(root)),sha256=sha(p)) for p in sorted(inputs)],
         outputs=[dict(path=str(p.relative_to(root)),sha256=sha(p)) for p in post.iterdir() if p.is_file()],
         processor_sha256=sha(Path(__file__)),failures=failures))
    try:accounting(root)
    except subprocess.CalledProcessError as error:(root/'accounting_error.txt').write_text(str(error))
    if complete:(root/'SUCCESS').write_text('All 18 cases and nine paired timing records validated.\n')
    else:(root/'FAILED').write_text('Incomplete or failed cases; see postprocessing/failures.json.\n')
    print(report/'REPORT.md')
    ensure(complete,'Campaign incomplete; valid partial data retained')


def archive(root):
    submitted=rows(root/'submitted_jobs.csv')
    target=[r['job_id'] for r in submitted if r['role']!='archive']
    result=subprocess.run(['squeue','-h','-u',getpass.getuser(),'-o','%i'],capture_output=True,text=True,check=True)
    ensure(not set(target).intersection(result.stdout.split()),'Compute/collector still active; refusing archive')
    result=subprocess.run(['sacct','-X','-n','-P','-j',','.join(target),'--format=JobIDRaw,State,ExitCode'],capture_output=True,text=True,check=True)
    found={row.split('|')[0]:row.split('|')[1:] for row in result.stdout.splitlines() if row.strip()}
    terminal={'COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED'}
    ensure(all(job in found and found[job][0].split()[0].rstrip('+') in terminal for job in target),'Accounting not terminal')
    accounting(root)
    dest=root.parent.parent/'results'/root.name
    ensure(not dest.exists(),'Archive destination already exists')
    save(root/'archive.json',dict(original_directory=str(root),archived_directory=str(dest),states=found,
                                 archive_job_id=os.environ.get('SLURM_JOB_ID'),archived_utc=datetime.now(timezone.utc).isoformat()))
    # The archive task never uses the campaign itself as working directory at rename.
    os.chdir(root.parent.parent)
    dest.parent.mkdir(exist_ok=True)
    root.rename(dest)
    if (dest/'source').is_dir() and read_manifest(dest).get('retain_source_after_completion') is False:
        ensure((dest/'source.sha256').is_file(), 'Source manifest missing')
        shutil.rmtree(dest/'source')
    print(dest,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['prepare','submit','build','case','collect','archive'])
    parser.add_argument('--run-root',type=Path)
    parser.add_argument('--dimension',choices=['1d','2d'])
    parser.add_argument('--restart',type=int,default=100)
    parser.add_argument('--minimum-memory-gib',type=int,default=16)
    parser.add_argument('--case')
    parser.add_argument('--cells',type=int,default=512)
    parser.add_argument('--node')
    args=parser.parse_args()
    if args.mode=='prepare':
        ensure(args.dimension is not None,'dimension required')
        print(prepare(args.dimension,args.restart,args.minimum_memory_gib,args.cells,args.node));return
    ensure(args.run_root is not None,'run-root required')
    root=args.run_root.resolve()
    if args.mode=='submit':submit(root)
    elif args.mode=='build':build(root)
    elif args.mode=='case':case(root,args.case)
    elif args.mode in ['collect','archive']:
        from provenance import runtime
        metadata=runtime(root/'bin/solver.exe' if (root/'bin/solver.exe').exists() else None)
        save(root/f'{args.mode}_environment.json',metadata)
        if args.mode=='collect':collect(root)
        else:archive(root)


if __name__=='__main__':main()
