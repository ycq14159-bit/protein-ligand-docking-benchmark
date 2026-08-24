#!/usr/bin/env bash
set -euo pipefail
ENV=/root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality/tool_envs/posebusters
LOG=/root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality/runs/20260812_full_01/logs/posebusters_install.log
if [[ ! -x "$ENV/bin/python" ]]; then
  /root/miniconda3/envs/entry-work-packages-v1/bin/python -m venv --system-site-packages "$ENV"
fi
"$ENV/bin/pip" install posebusters >"$LOG" 2>&1
"$ENV/bin/python" - <<'PY'
import importlib.metadata
import posebusters
print(importlib.metadata.version("posebusters"))
PY
