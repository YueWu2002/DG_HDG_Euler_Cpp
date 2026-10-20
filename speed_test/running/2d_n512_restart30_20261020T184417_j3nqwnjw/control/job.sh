#!/usr/bin/env bash
set -euo pipefail
module purge
module load intel-oneapi-compilers intel-oneapi-mkl eigen cmake
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
exec srun --export=ALL --ntasks=1 --cpus-per-task=1 --cpu-bind=cores python3 "$@"
