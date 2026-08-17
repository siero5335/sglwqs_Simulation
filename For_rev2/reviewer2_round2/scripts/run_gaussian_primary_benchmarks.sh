#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

RSCRIPT="${RSCRIPT:-Rscript}"
export R2R2_REPO_ROOT="${R2R2_REPO_ROOT:-$(pwd)}"
export R2R2_ROOT="${R2R2_ROOT:-$(pwd)/reviewer2_round2}"
export R2R2_SGLWQS_SOURCE="${R2R2_SGLWQS_SOURCE:-$(pwd)/source/sglwqs-envint-revision-docs}"
export R2R2_WORKERS="${R2R2_WORKERS:-15}"
export R2R2_N_BOOT="${R2R2_N_BOOT:-200}"

"$RSCRIPT" reviewer2_round2/scripts/60_gaussian_primary_benchmarks.R
