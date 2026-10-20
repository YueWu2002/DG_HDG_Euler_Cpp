#!/usr/bin/env python3
"""Move only scheduler-confirmed terminal runs from running/ to results/."""
import argparse
import json
from pathlib import Path
import subprocess
import shutil
import getpass
from datetime import datetime, timezone

TERMINAL = {'COMPLETED', 'FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY',
            'NODE_FAIL', 'PREEMPTED', 'BOOT_FAIL', 'DEADLINE', 'REVOKED'}


def archive(root, only=None):
    for folder in sorted((root/'running').iterdir()):
        if only is not None and folder.resolve() != only.resolve():
            continue
        if not folder.is_dir() or not (folder/'job_id.txt').exists():
            continue
        job = (folder/'job_id.txt').read_text().strip().split(';')[0]
        if not job.isdigit():
            raise ValueError(f'Invalid job ID: {folder}')
        live = subprocess.run(['squeue', '-h', '-j', job, '-o', '%i'], capture_output=True, text=True)
        if live.returncode or live.stdout.strip():
            continue
        accounting = subprocess.run(['sacct', '-X', '-j', job, '-n', '-P',
                                     '--format=JobIDRaw,State,ExitCode'], check=True, capture_output=True, text=True)
        matches = [line.split('|') for line in accounting.stdout.splitlines() if line.split('|')[0] == job]
        if len(matches) != 1 or matches[0][1].split()[0].rstrip('+') not in TERMINAL:
            continue
        # Other active jobs may still use this run's scripts or snapshot.
        active = subprocess.run(['squeue', '-h', '-u', getpass.getuser(), '-o', '%E|%Z'],
                                capture_output=True, text=True)
        if active.returncode:
            continue
        if any(str(folder.resolve()) in line or job in line.split('|')[0].replace(':', ',').split(',')
               for line in active.stdout.splitlines()):
            continue
        destination = root/'results'/folder.name
        if destination.exists():
            raise FileExistsError(destination)
        record = dict(job_id=job, state=matches[0][1], exit_code=matches[0][2],
                      archived_utc=datetime.now(timezone.utc).isoformat(), original_directory=str(folder.resolve()))
        (folder/'archive.json').write_text(json.dumps(record, indent=2)+'\n')
        destination.parent.mkdir(parents=True, exist_ok=True)
        folder.rename(destination)
        manifest = destination/'submission.json'
        if manifest.exists() and json.loads(manifest.read_text()).get('retain_source_after_completion') is False:
            source = destination/'source'
            if source.is_symlink():
                raise RuntimeError(f'Refusing source symlink: {source}')
            if source.is_dir():
                if not (destination/'source.sha256').is_file():
                    raise RuntimeError('Source manifest required before removing temporary source')
                shutil.rmtree(source)
                record['source_removed'] = True
                record['source_removal_policy'] = 'Keep SHA256 manifest; remove temporary source after completion'
                (destination/'archive.json').write_text(json.dumps(record, indent=2)+'\n')

        print(f'{folder} -> {destination}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--run-root', type=Path, help='Archive only this scheduler-confirmed finished run')
    args = parser.parse_args()
    for name in ('accuracy_test', 'speed_test'):
        root = args.root/name
        (root/'running').mkdir(parents=True, exist_ok=True)
        (root/'results').mkdir(parents=True, exist_ok=True)
        archive(root, args.run_root)
