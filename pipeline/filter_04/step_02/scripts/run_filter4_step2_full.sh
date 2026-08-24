#!/usr/bin/env bash
set -euo pipefail

cd /root/autodl-tmp/benchmark_1.0
exec /root/miniconda3/envs/entry-work-packages-v1/bin/python \
  filter_04_crystal_packing_influence/step_02_biological_assembly_equivalence/scripts/filter4_step2_pipeline.py \
  --config filter_04_crystal_packing_influence/step_02_biological_assembly_equivalence/config/filter4_step2_v1.yaml \
  --run-dir filter_04_crystal_packing_influence/step_02_biological_assembly_equivalence/runs/step02_full_v3 \
  --mode full
