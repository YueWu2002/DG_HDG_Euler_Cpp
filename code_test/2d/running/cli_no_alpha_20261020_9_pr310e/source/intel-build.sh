#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
module purge
module load intel-oneapi-compilers intel-oneapi-mkl eigen cmake
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1
cmake -S "$PROJECT_ROOT" -B "${BUILD_DIR:-$PROJECT_ROOT/build-intel}" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=icpx \
  -DDGH_NATIVE_ARCH=ON -DDGH_BLAS_BACKEND=MKL -DDGH_LAPACK_BACKEND=MKL \
  -DDGH_BUILD_1D="${BUILD_1D:-ON}" -DDGH_BUILD_2D="${BUILD_2D:-OFF}" \
  -DBUILD_TESTING="${BUILD_TESTING:-ON}"
cmake --build "${BUILD_DIR:-$PROJECT_ROOT/build-intel}" -j "${BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-4}}"
