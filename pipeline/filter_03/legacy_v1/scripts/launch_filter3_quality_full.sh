#!/usr/bin/env bash
set -euo pipefail
source /root/miniconda3/etc/profile.d/conda.sh
conda activate entry-work-packages-v1
exec python /root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality/scripts/filter3_quality_pipeline.py \
  --start 0 --end 255 --workers 4
