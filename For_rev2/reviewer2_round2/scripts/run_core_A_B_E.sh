#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

RSCRIPT="${RSCRIPT:-Rscript}"
export R2R2_REPO_ROOT="${R2R2_REPO_ROOT:-$(pwd)}"
export R2R2_ROOT="${R2R2_ROOT:-$(pwd)/reviewer2_round2}"
export R2R2_WORKERS="${R2R2_WORKERS:-8}"
export R2R2_N_BOOT="${R2R2_N_BOOT:-200}"
LEGACY_SOURCE="${R2R2_SGLWQS_LEGACY_SOURCE:-$(pwd)/source/sglwqs_old_analysis/sglwqs-main}"
REVISION_SOURCE="${R2R2_SGLWQS_REVISION_SOURCE:-$(pwd)/source/sglwqs-envint-revision-docs}"

R2R2_SGLWQS_SOURCE="$REVISION_SOURCE" "$RSCRIPT" reviewer2_round2/scripts/00_audit_existing_simulations.R
R2R2_SGLWQS_SOURCE="$REVISION_SOURCE" "$RSCRIPT" reviewer2_round2/scripts/01_smoke_test.R
R2R2_SGLWQS_SOURCE="$LEGACY_SOURCE" "$RSCRIPT" reviewer2_round2/scripts/10_gaussian_sample_size.R
R2R2_SGLWQS_SOURCE="$REVISION_SOURCE" "$RSCRIPT" reviewer2_round2/scripts/20_gaussian_high_dimensional.R
R2R2_SGLWQS_SOURCE="$REVISION_SOURCE" "$RSCRIPT" reviewer2_round2/scripts/50_hard_setting_sample_size.R
