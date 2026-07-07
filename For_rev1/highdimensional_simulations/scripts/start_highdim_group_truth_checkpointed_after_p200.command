#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}" || exit 1
mkdir -p logs

WAIT_PID="${1:-${GT_WAIT_PID:-14047}}"

bash analysis/wait_then_start_highdim_group_truth_checkpointed.sh "${WAIT_PID}" 2>&1 | tee logs/highdim_group_truth_checkpointed_after_p200.log
