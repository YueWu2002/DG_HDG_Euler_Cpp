#!/usr/bin/env bash
#SBATCH --job-name=eigen-euler-2d
#SBATCH --partition=batch
#SBATCH --constraint=intel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=1-00:00:00
set -euo pipefail
: "${DGH_EXPERIMENT:?}" "${DGH_RESULT_ROOT:?}" "${SLURM_JOB_ID:?}"
RUN_ROOT="$DGH_RESULT_ROOT/job_$SLURM_JOB_ID"
mkdir "$RUN_ROOT"
python3 -c 'import sys; sys.path.insert(0,"scripts"); from provenance import runtime,save; save(sys.argv[1],runtime())' "$RUN_ROOT/provenance.json"
trap 'status=$?; printf "%s\n" "$status" > "$RUN_ROOT/exit_code.txt"' EXIT
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1
export BUILD_DIR="$RUN_ROOT/build" BUILD_JOBS=1 BUILD_TESTING=ON BUILD_1D=OFF BUILD_2D=ON
bash intel-build.sh > "$RUN_ROOT/build.log" 2>&1
module purge
module load intel-oneapi-compilers intel-oneapi-mkl eigen cmake
PYTHON_EXE=$(command -v python3)
{
  date -Is
  hostname
  lscpu
  icpx --version
  cmake --version
  module list
  sha256sum "$BUILD_DIR/dg_hdg_2d/dg_hdg_2d"
  printf 'python=%s\n' "$PYTHON_EXE"
} > "$RUN_ROOT/environment.txt" 2>&1
ctest --test-dir "$BUILD_DIR" --output-on-failure > "$RUN_ROOT/ctest.log" 2>&1
"$PYTHON_EXE" code_test/2d/test_solver.py --exe "$BUILD_DIR/dg_hdg_2d/dg_hdg_2d" --output "$RUN_ROOT/checks" > "$RUN_ROOT/code_test.log" 2>&1
if [[ "$DGH_EXPERIMENT" == accuracy ]]; then
  srun --export=ALL --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
    "$PYTHON_EXE" accuracy_test/2d/convergence.py --exe "$BUILD_DIR/dg_hdg_2d/dg_hdg_2d" --output "$RUN_ROOT/data"
else
  srun --export=ALL --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
    "$PYTHON_EXE" speed_test/2d/benchmark.py --exe "$BUILD_DIR/dg_hdg_2d/dg_hdg_2d" --output "$RUN_ROOT/data"
fi
