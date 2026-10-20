#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
exec srun --ntasks=1 --cpus-per-task=1 /usr/bin/python3 "$@"
