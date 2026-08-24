#!/usr/bin/env bash
set -euo pipefail
source /root/miniconda3/etc/profile.d/conda.sh
conda activate entry-work-packages-v1
ROOT=/root/autodl-tmp/benchmark_1.0/filter_04_crystal_packing_influence/step_01_lattice_neighbor_search
mkdir -p "$ROOT/runs/step01_full_v3"
PILOT="$ROOT/runs/step01_pilot_v3/validation.json"
python - "$PILOT" <<'PY'
import json, sys
value = json.load(open(sys.argv[1]))
if not value.get('validation_pass'):
    raise SystemExit('pilot validation is not PASS')
PY
python "$ROOT/scripts/filter4_step1_pipeline.py" \
  --config "$ROOT/config/filter4_step1.yaml" \
  --run-dir "$ROOT/runs/step01_full_v3" \
  --mode full --workers 7 \
  > "$ROOT/runs/step01_full_v3/runtime.log" 2>&1
