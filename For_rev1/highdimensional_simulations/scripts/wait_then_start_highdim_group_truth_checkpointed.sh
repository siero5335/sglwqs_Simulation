#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
mkdir -p logs

WAIT_PID="${1:-${GT_WAIT_PID:-14047}}"
SLEEP_SECONDS="${GT_WAIT_SLEEP_SECONDS:-300}"

echo "[waiter] started at $(date)"
echo "[waiter] waiting for PID ${WAIT_PID}; polling every ${SLEEP_SECONDS}s"

while kill -0 "${WAIT_PID}" 2>/dev/null; do
  echo "[waiter] $(date): PID ${WAIT_PID} is still running"
  sleep "${SLEEP_SECONDS}"
done

echo "[waiter] $(date): PID ${WAIT_PID} is no longer running"
echo "[waiter] launching highdim group/truth checkpointed run"

bash ./run_highdim_group_truth_checkpointed.sh

echo "[waiter] finished at $(date)"
