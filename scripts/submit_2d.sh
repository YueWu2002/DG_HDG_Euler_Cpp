#!/usr/bin/env bash
set -euo pipefail
EXPERIMENT=${1:?accuracy or speed}
shift
[[ "$EXPERIMENT" == accuracy || "$EXPERIMENT" == speed ]]
PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
RESULT_ROOT="$PROJECT_ROOT/${EXPERIMENT}_test/running"
mkdir -p "$RESULT_ROOT" "$PROJECT_ROOT/${EXPERIMENT}_test/results"
SUBMISSION=$(mktemp -d "$RESULT_ROOT/2d_submission_$(date +%Y%m%dT%H%M%S)_XXXXXX")
SNAPSHOT="$SUBMISSION/source"
mkdir "$SNAPSHOT"
tar -C "$PROJECT_ROOT" --exclude=__pycache__ --exclude=results --exclude=running -cf - CMakeLists.txt cmake common dg_hdg_1d dg_hdg_2d vendor scripts \
  intel-build.sh skills AGENTS.md README.md MATHEMATICS.md DESIGN.md VALIDATION.md MEMORY_ESTIMATE.md RESULTS.md RESULTS_2D.md \
  code_test/2d accuracy_test/2d speed_test/2d | tar -C "$SNAPSHOT" -xf -
git -C "$PROJECT_ROOT" rev-parse HEAD > "$SUBMISSION/base_commit.txt"
git -C "$PROJECT_ROOT" diff --binary > "$SUBMISSION/working_tree.patch"
(cd "$SNAPSHOT" && find . -type f -print0 | sort -z | xargs -0 sha256sum) > "$SUBMISSION/source.sha256"
python3 "$PROJECT_ROOT/scripts/provenance.py" "$SUBMISSION" 2d "$EXPERIMENT" -- "$@"
MEMORY=100G
if [[ "$EXPERIMENT" == accuracy ]]; then MEMORY=16G; fi
JOB_ID=$(sbatch --mem="$MEMORY" --parsable --chdir="$SNAPSHOT" --job-name="eigen-2d-$EXPERIMENT" \
  --export="ALL,DGH_EXPERIMENT=$EXPERIMENT,DGH_RESULT_ROOT=$SUBMISSION" \
  --output="$SUBMISSION/%x-%j.out" --error="$SUBMISSION/%x-%j.err" "$@" "$SNAPSHOT/scripts/experiment_2d_job.sh")
printf '%s\n' "$JOB_ID" > "$SUBMISSION/job_id.txt"
printf 'job_id=%s\noutput=%s\n' "$JOB_ID" "$SUBMISSION"
