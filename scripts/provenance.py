"""Persistent source and runtime identity for experiment artifacts (stdlib only)."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
from datetime import datetime, timezone


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def save(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
    temporary.replace(path)


def source_identity(root):
    root = Path(root).resolve()
    for parent in (root, *root.parents):
        manifest = parent / 'submission.json'
        if manifest.is_file():
            return json.loads(manifest.read_text())
    # Local runs have no submission snapshot: hash all solver/build sources.
    files = []
    for name in ('CMakeLists.txt', 'cmake', 'common', 'dg_hdg_1d', 'dg_hdg_2d', 'vendor', 'scripts', 'intel-build.sh'):
        path = root / name
        files.extend([path] if path.is_file() else path.rglob('*') if path.exists() else [])
    h = hashlib.sha256()
    for path in sorted(p for p in files if p.is_file() and '__pycache__' not in p.parts):
        h.update(str(path.relative_to(root)).encode() + b'\0' + digest(path).encode() + b'\n')
    commit = subprocess.run(['git', '-C', str(root), 'rev-parse', 'HEAD'], capture_output=True, text=True)
    sha = h.hexdigest()
    base = commit.stdout.strip() or 'unknown'
    return dict(code_version=f'{base[:12]}+src.{sha[:16]}', base_commit=base, source_sha256=sha)


def runtime(exe=None):
    result = source_identity(Path(__file__).resolve().parents[1])
    cpu = {}
    if Path('/proc/cpuinfo').exists():
        for line in Path('/proc/cpuinfo').read_text().splitlines():
            if ':' in line:
                key, value = line.split(':', 1)
                if key.strip() in ('model name', 'vendor_id', 'cpu family', 'model', 'stepping', 'microcode'):
                    cpu.setdefault(key.strip(), value.strip())
    result.update(recorded_utc=datetime.now(timezone.utc).isoformat(), host=platform.node(),
                  cpu=cpu, cpu_model=cpu.get('model name', platform.processor()),
                  cpu_affinity=sorted(os.sched_getaffinity(0)),
                  slurm={k: v for k, v in os.environ.items() if k.startswith('SLURM_')},
                  threads={k: os.environ.get(k, '') for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'BLIS_NUM_THREADS')})
    if exe:
        exe = Path(exe).resolve()
        result.update(executable=str(exe), executable_sha256=digest(exe))
        for parent in exe.parents:
            cache = parent / 'CMakeCache.txt'
            if cache.exists():
                result['cmake_cache'] = cache.read_text()
                break
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('submission', type=Path)
    parser.add_argument('dimension')
    parser.add_argument('experiment')
    parser.add_argument('sbatch_args', nargs='*')
    args = parser.parse_args()
    sha = digest(args.submission / 'source.sha256')
    base = (args.submission / 'base_commit.txt').read_text().strip()
    save(args.submission / 'submission.json', dict(schema_version=1,
         code_version=f'{base[:12]}+src.{sha[:16]}', base_commit=base, source_sha256=sha,
         dimension=args.dimension, experiment=args.experiment,
         submitted_utc=datetime.now(timezone.utc).isoformat(), sbatch_args=args.sbatch_args,
         source_manifest='source.sha256', source_directory='source',
         retain_source_after_completion=False,
         parameters='job_<id>/data/*.command.json', runtime='job_<id>/data/provenance.json'))
