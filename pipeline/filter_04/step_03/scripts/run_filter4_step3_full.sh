#!/usr/bin/env bash
set -euo pipefail

cd /root/autodl-tmp/benchmark_1.0
exec /root/miniconda3/envs/entry-work-packages-v1/bin/python \
  filter_04_crystal_packing_influence/step_03_direct_ligand_crystal_contact/scripts/filter4_step3_pipeline.py \
  --config filter_04_crystal_packing_influence/step_03_direct_ligand_crystal_contact/config/filter4_step3_v1.yaml \
  --run-dir filter_04_crystal_packing_influence/step_03_direct_ligand_crystal_contact/runs/step03_full_v2 \
  --mode full
