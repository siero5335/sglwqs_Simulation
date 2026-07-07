#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}" || exit 1
mkdir -p logs

bash ./run_highdim_p_checkpointed_10core.sh 2>&1 | tee logs/highdim_p_checkpointed_10core.log
