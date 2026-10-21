#!/usr/bin/env bash
# Called from a source snapshot on a compute node before a Release experiment.
set -euo pipefail
dimension=${1:?1d or 2d}
run_root=${2:?run directory}
[[ "$dimension" == 1d || "$dimension" == 2d ]]
project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export BUILD_DIR="$run_root/build-debug" BUILD_TYPE=Debug BUILD_TESTING=ON BUILD_JOBS=1
export BUILD_1D=OFF BUILD_2D=OFF
if [[ "$dimension" == 1d ]]; then export BUILD_1D=ON; solver=dg_hdg; else export BUILD_2D=ON; solver=dg_hdg_2d; fi
bash "$project_root/intel-build.sh" > "$run_root/build-debug.log" 2>&1
module purge
module load intel-oneapi-compilers intel-oneapi-mkl eigen cmake
ctest --test-dir "$BUILD_DIR" --output-on-failure > "$run_root/ctest.log" 2>&1
DGH_CHECK_ROOT="$run_root" python3 "$project_root/code_test/$dimension/test_solver.py" --exe "$BUILD_DIR/bin/$solver.exe" \
  --output "$run_root/checks" > "$run_root/code_test.log" 2>&1
