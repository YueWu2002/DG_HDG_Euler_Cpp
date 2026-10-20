#!/usr/bin/env bash
#SBATCH --partition=batch
#SBATCH --constraint=intel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
set -euo pipefail
RUN_ROOT=$(dirname "$PWD")
trap 'rc=$?; printf "%s\n" "$rc" > "$RUN_ROOT/exit_code.txt"' EXIT
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1
python3 -c 'import sys;sys.path.insert(0,"scripts");from provenance import runtime,save;save(sys.argv[1],runtime())' "$RUN_ROOT/provenance.json"
export BUILD_DIR="$RUN_ROOT/build" BUILD_JOBS=2 BUILD_1D=OFF BUILD_2D=ON
bash intel-build.sh > "$RUN_ROOT/build.log" 2>&1
module purge
module load intel-oneapi-compilers intel-oneapi-mkl eigen cmake
{ hostname; lscpu; icpx --version; module list; } > "$RUN_ROOT/environment.txt" 2>&1
ctest --test-dir "$BUILD_DIR" --output-on-failure > "$RUN_ROOT/ctest.log" 2>&1
python3 code_test/2d/test_solver.py --exe "$BUILD_DIR/dg_hdg_2d/dg_hdg_2d" --output "$RUN_ROOT/checks" > "$RUN_ROOT/code_test.log" 2>&1
python3 - "$BUILD_DIR/dg_hdg_2d/dg_hdg_2d" "$RUN_ROOT" <<'CHECK'
import subprocess,sys,json
from pathlib import Path
exe=sys.argv[1];root=Path(sys.argv[2])
help_result=subprocess.run([exe,'--help'],capture_output=True,text=True,check=True)
(root/'help.txt').write_text(help_result.stdout)
assert 'alpha' not in help_result.stdout
rejected=subprocess.run([exe,'--alpha-deg','30'],capture_output=True,text=True)
(root/'removed_option.log').write_text(rejected.stdout+rejected.stderr)
assert rejected.returncode!=0 and 'Unknown option: --alpha-deg' in rejected.stderr
(root/'SUCCESS').write_text('CTest, production checks, help and removed option checks passed.\n')
CHECK
