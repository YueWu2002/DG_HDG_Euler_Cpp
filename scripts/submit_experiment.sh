#!/usr/bin/env bash
set -euo pipefail
EXPERIMENT=${1:?accuracy or speed}
shift
[[ "$EXPERIMENT" == accuracy || "$EXPERIMENT" == speed ]]
PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
RESULT_ROOT="$PROJECT_ROOT/${EXPERIMENT}_test/running"
mkdir -p "$RESULT_ROOT" "$PROJECT_ROOT/${EXPERIMENT}_test/results"
SUBMISSION=$(mktemp -d "$RESULT_ROOT/1d_submission_$(date +%Y%m%dT%H%M%S)_XXXXXX")
SNAPSHOT="$SUBMISSION/source"
mkdir "$SNAPSHOT"
tar -C "$PROJECT_ROOT" --exclude=__pycache__ --exclude=results --exclude=running -cf - CMakeLists.txt cmake common dg_hdg_1d vendor scripts \
  intel-build.sh skills AGENTS.md README.md MATHEMATICS.md DESIGN.md VALIDATION.md MEMORY_ESTIMATE.md RESULTS.md \
  code_test/1d accuracy_test/1d speed_test/1d | tar -C "$SNAPSHOT" -xf -
git -C "$PROJECT_ROOT" rev-parse HEAD > "$SUBMISSION/base_commit.txt"
git -C "$PROJECT_ROOT" diff --binary > "$SUBMISSION/working_tree.patch"
(cd "$SNAPSHOT" && find . -type f -print0 | sort -z | xargs -0 sha256sum) > "$SUBMISSION/source.sha256"
python3 "$PROJECT_ROOT/scripts/provenance.py" "$SUBMISSION" 1d "$EXPERIMENT" -- "$@"
JOB_ID=$(sbatch --parsable --chdir="$SNAPSHOT" --job-name="eigen-1d-$EXPERIMENT" \
  --export="ALL,DGH_EXPERIMENT=$EXPERIMENT,DGH_RESULT_ROOT=$SUBMISSION" \
  --output="$SUBMISSION/%x-%j.out" --error="$SUBMISSION/%x-%j.err" \
  "$@" "$SNAPSHOT/scripts/experiment_job.sh")
printf '%s\n' "$JOB_ID" > "$SUBMISSION/job_id.txt"
printf 'job_id=%s\noutput=%s\n' "$JOB_ID" "$SUBMISSION"
