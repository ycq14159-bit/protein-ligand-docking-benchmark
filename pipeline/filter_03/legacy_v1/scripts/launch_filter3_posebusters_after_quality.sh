#!/usr/bin/env bash
set -euo pipefail
RUN=/root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality/runs/20260812_full_01
while true; do
  completed=$(find "$RUN/work/quality_batches" -mindepth 2 -maxdepth 2 -name _COMPLETE.json | wc -l)
  if [[ "$completed" -ge 256 ]]; then
    break
  fi
  printf '%s waiting quality buckets: %s/256\n' "$(date -Is)" "$completed"
  sleep 60
done
exec /root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality/tool_envs/posebusters/bin/python \
  /root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality/scripts/filter3_posebusters.py \
  --start 0 --end 255 --workers 8
