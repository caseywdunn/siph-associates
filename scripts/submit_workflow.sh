#!/usr/bin/env bash
set -euo pipefail

PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
TARGET="${1:-phase1_ready}"
CONFIG="${2:-config/config.yaml}"
mkdir -p "$PROJECT/logs/slurm" "$PROJECT/workflow/state"

JOB_ID=$(sbatch --parsable \
    --export="ALL,SIPH_ASSOCIATES_PROJECT=$PROJECT,SIPH_ASSOCIATES_TARGET=$TARGET,SIPH_ASSOCIATES_CONFIG=$CONFIG" \
    "$PROJECT/scripts/workflow_controller.sbatch")
JOB_ID="${JOB_ID%%;*}"
COMMIT=$(git -C "$PROJECT" rev-parse HEAD)
TEMP="$PROJECT/workflow/state/latest_controller.tsv.tmp"
printf 'submitted_at\tcontroller_job_id\ttarget\tconfig\tcommit\n' > "$TEMP"
printf '%s\t%s\t%s\t%s\t%s\n' "$(date -Iseconds)" "$JOB_ID" "$TARGET" "$CONFIG" "$COMMIT" >> "$TEMP"
mv "$TEMP" "$PROJECT/workflow/state/latest_controller.tsv"
printf '%s\n' "$JOB_ID"
