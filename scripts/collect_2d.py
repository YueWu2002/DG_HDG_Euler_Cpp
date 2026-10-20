#!/usr/bin/env python3
"""Compatibility entry: postprocess each experiment in its own archive."""
import argparse
from pathlib import Path
import sys
from postprocess import accuracy, speed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--accuracy-job-root', type=Path, required=True)
    parser.add_argument('--speed-job-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, help='Deprecated: reports now belong to their respective run roots')
    args = parser.parse_args()
    if args.output:
        print('--output is deprecated; no combined data/report copy will be created.', file=sys.stderr)
    for job_root, processor in [(args.accuracy_job_root, accuracy), (args.speed_job_root, speed)]:
        root = job_root.resolve()
        if not root.name.startswith('job_'):
            parser.error('Expected a job_<id> directory')
        processor(root.parent)


if __name__ == '__main__':
    main()
