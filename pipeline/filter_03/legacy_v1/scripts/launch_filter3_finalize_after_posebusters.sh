#!/usr/bin/env bash
set -euo pipefail
RUN=/root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality/runs/20260812_full_01
while true; do
  completed=$(find "$RUN/work/posebusters_batches" -mindepth 2 -maxdepth 2 -name _COMPLETE.json | wc -l)
  if [[ "$completed" -ge 256 ]]; then
    break
  fi
  printf '%s waiting PoseBusters buckets: %s/256\n' "$(date -Is)" "$completed"
  sleep 60
done
source /root/miniconda3/etc/profile.d/conda.sh
conda activate entry-work-packages-v1
exec python /root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality/scripts/filter3_finalize.py
