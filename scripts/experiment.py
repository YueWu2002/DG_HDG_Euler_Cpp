"""Small shared helpers for reproducible, single-thread 1D experiments."""
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import subprocess
import time
import tempfile
from datetime import datetime
from provenance import runtime, save

_CONTEXT = {}


def rows(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, data):
    if not data:
        return
    data = [dict(row, **_CONTEXT) for row in data]
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)


def run(command, output):
    command = [str(part) for part in command]
    output = Path(output)
    output.with_suffix('.command.json').write_text(json.dumps(command, indent=2) + '\n')
    env = dict(os.environ)
    for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        env[key] = '1'
    metadata = runtime(command[0])
    metadata['command'] = command
    metadata['parameters'] = command[1:]
    save(output.with_suffix('.provenance.json'), metadata)
    start = time.perf_counter()
    # wait4 records this child's CPU time and RSS, rather than a cumulative
    # maximum across previous DG/HDG children.
    with output.with_suffix('.log').open('w') as stream:
        process = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT, env=env)
        _, status, usage = os.wait4(process.pid, 0)
        process.returncode = os.waitstatus_to_exitcode(status)
    metrics = dict(returncode=process.returncode, wall_s=time.perf_counter()-start,
                   user_cpu_s=usage.ru_utime, system_cpu_s=usage.ru_stime,
                   peak_rss_kib=usage.ru_maxrss, host=platform.node(),
                   slurm_job_id=os.environ.get('SLURM_JOB_ID', ''))
    output.with_suffix('.resources.json').write_text(json.dumps(metrics, indent=2) + '\n')
    if process.returncode:
        raise RuntimeError(f'Command failed: {command}; see {output.with_suffix(".log")}')
    return metrics


def require_debug_build(exe):
    """code_test requires a separately configured Debug solver."""
    exe = Path(exe).resolve(strict=True)
    for parent in exe.parents:
        cache = parent/'CMakeCache.txt'
        if cache.is_file():
            if 'CMAKE_BUILD_TYPE:STRING=Debug' not in cache.read_text().splitlines():
                raise RuntimeError('code_test requires a Debug build; use BUILD_TYPE=Debug with a separate BUILD_DIR')
            return exe
    raise RuntimeError('Cannot verify Debug build: CMakeCache.txt not found for code_test executable')


def prepare(output, exe):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise RuntimeError(f'Output must be empty: {output}')
    exe = Path(exe).resolve(strict=True)
    (output/'executable.json').write_text(json.dumps(dict(
        path=str(exe), sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
        host=platform.node(), slurm_job_id=os.environ.get('SLURM_JOB_ID', '')), indent=2)+'\n')
    metadata = runtime(exe)
    save(output/'provenance.json', metadata)
    _CONTEXT.clear()
    _CONTEXT.update(code_version=metadata['code_version'], source_sha256=metadata['source_sha256'],
                    executable_sha256=metadata['executable_sha256'], cpu_model=metadata['cpu_model'],
                    host=metadata['host'], slurm_job_id=os.environ.get('SLURM_JOB_ID', ''))
    return output, exe


def prepare_experiment(output, exe, experiment, dimension):
    """Keep local and batch experiment outputs inside the live run tree."""
    project = Path(__file__).resolve().parents[1]
    submission = os.environ.get('DGH_RESULT_ROOT')
    running = (Path(submission).resolve().parent if submission else
               project/f'{experiment}_test'/'running')
    if running.name != 'running':
        raise RuntimeError(f'New experiments must run below running/: {running}')
    # Validate executable before creating a new run directory.
    exe = Path(exe).resolve(strict=True)
    if output is None:
        running.mkdir(parents=True, exist_ok=True)
        output = Path(tempfile.mkdtemp(
            prefix=f'{dimension}_local_{datetime.now():%Y%m%dT%H%M%S}_', dir=running))
    else:
        output = Path(output).resolve()
        if running not in output.parents:
            raise RuntimeError(f'--output must be inside {running}; got {output}')
    print(f'Output: {output}', flush=True)
    return prepare(output, exe)
