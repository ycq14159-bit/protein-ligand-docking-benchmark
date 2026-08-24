#!/usr/bin/env bash
set -euo pipefail
cd /root/autodl-tmp/benchmark_1.0/filter_2_ligand_qualification_v3
exec /root/miniconda3/envs/interaction-pilot-v2/bin/python \
  scripts/filter2_v3_pipeline.py run-full --workers 24 --batch-size 200 \
  >> runs/20260804_full_01/logs/full.screen.log 2>&1
