#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality
RUN="$ROOT/runs/20260812_full_01"
exec /root/miniconda3/envs/entry-work-packages-v1/bin/python \
  "$ROOT/scripts/filter3_pipeline.py" download --workers 24 \
  >> "$RUN/logs/validation_download.log" 2>&1
