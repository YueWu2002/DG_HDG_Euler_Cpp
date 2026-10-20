#!/usr/bin/env bash
set -euo pipefail
ACCURACY_ROOT=$(realpath -m "${1:?accuracy job directory}")
SPEED_ROOT=$(realpath -m "${2:?speed job directory}")
ACCURACY_ID=$(cat "$(dirname "$ACCURACY_ROOT")/job_id.txt")
SPEED_ID=$(cat "$(dirname "$SPEED_ROOT")/job_id.txt")
[[ "$ACCURACY_ID" =~ ^[0-9]+$ && "$SPEED_ID" =~ ^[0-9]+$ ]]
[[ "$(basename "$ACCURACY_ROOT")" == "job_$ACCURACY_ID" && "$(basename "$SPEED_ROOT")" == "job_$SPEED_ID" ]]
PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
mkdir -p "$PROJECT_ROOT/speed_test/running"
REPORT_ROOT=$(mktemp -d "$PROJECT_ROOT/speed_test/running/2d_report_submission_$(date +%Y%m%dT%H%M%S)_XXXXXX")
cp "$PROJECT_ROOT/scripts/collect_2d.py" "$PROJECT_ROOT/scripts/postprocess.py" "$REPORT_ROOT/"
sha256sum "$REPORT_ROOT/collect_2d.py" "$REPORT_ROOT/postprocess.py" > "$REPORT_ROOT/collector.sha256"
{
    printf '#!/usr/bin/env bash\nset -euo pipefail\n'
    printf '/usr/bin/python3 %q --accuracy-job-root %q --speed-job-root %q\n' \
        "$REPORT_ROOT/collect_2d.py" "$ACCURACY_ROOT" "$SPEED_ROOT"
    printf 'sacct -j %q --format=JobID,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList > %q\n' \
        "$ACCURACY_ID,$SPEED_ID" "$REPORT_ROOT/slurm_accounting.txt"
} > "$REPORT_ROOT/run.sh"
JOB_ID=$(sbatch --parsable --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=1G --time=00:05:00 \
    --job-name=eigen-2d-report --dependency="afterok:$ACCURACY_ID:$SPEED_ID" --chdir="$REPORT_ROOT" \
    --output="$REPORT_ROOT/%x-%j.out" --error="$REPORT_ROOT/%x-%j.err" "$REPORT_ROOT/run.sh")
printf '%s\n' "$JOB_ID" > "$REPORT_ROOT/job_id.txt"
printf 'job_id=%s\npostprocessing_log=%s\n' "$JOB_ID" "$REPORT_ROOT"
