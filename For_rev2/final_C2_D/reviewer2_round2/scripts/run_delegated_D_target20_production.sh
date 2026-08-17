#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ARCHIVE_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
WORK_ROOT="$ARCHIVE_ROOT/final_C2_D"
RSCRIPT="${RSCRIPT:-Rscript}"

cd "$WORK_ROOT"
unset R2R2_RESULT_NAMESPACE
export R2R2_REPO_ROOT="$WORK_ROOT"
export R2R2_ROOT="$WORK_ROOT/reviewer2_round2"
export R2R2_SGLWQS_SOURCE="$ARCHIVE_ROOT/source/sglwqs-envint-revision-docs"
export R2R2_SGLWQS_SOURCE_COMMIT=2fdd519e520a7dad1162810643e175cd616b1154
export R2R2_WORKERS="${R2R2_WORKERS:-8}"
export R2R2_N_BOOT="${R2R2_N_BOOT:-200}"
export R2R2_FORCE="${R2R2_FORCE:-false}"
export R2R2_SAVE_RAW_FIT="${R2R2_SAVE_RAW_FIT:-false}"

"$RSCRIPT" reviewer2_round2/scripts/04_preflight_delegated_D_target20.R
"$RSCRIPT" reviewer2_round2/scripts/40_unbalanced_effect_factorial.R
R2R2_VALIDATE_BATTERY=D "$RSCRIPT" reviewer2_round2/scripts/03_validate_delegated_battery.R
"$RSCRIPT" reviewer2_round2/scripts/41_summarize_D_target20.R
