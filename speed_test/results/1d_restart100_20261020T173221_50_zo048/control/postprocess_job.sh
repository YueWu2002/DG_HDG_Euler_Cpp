#!/usr/bin/env bash
set -euo pipefail
exec srun --ntasks=1 --cpus-per-task=1 /usr/bin/python3 "$@"
