#!/usr/bin/env bash
set -euo pipefail
source /root/miniconda3/etc/profile.d/conda.sh
conda activate entry-work-packages-v1
ROOT=/root/autodl-tmp/benchmark_1.0/filter_04_crystal_packing_influence/step_01_lattice_neighbor_search
mkdir -p "$ROOT/runs/step01_pilot_v3"
python "$ROOT/scripts/filter4_step1_pipeline.py" \
  --config "$ROOT/config/filter4_step1.yaml" \
  --run-dir "$ROOT/runs/step01_pilot_v3" \
  --mode pilot --workers 3 \
  > "$ROOT/runs/step01_pilot_v3/runtime.log" 2>&1
