#!/usr/bin/env bash
set -euo pipefail

PY=/root/miniconda3/envs/interaction-pilot-v2/bin/python
SCRIPT=/root/autodl-tmp/benchmark_1.0/filter_2_ligand_qualification_v2/scripts/filter2_v2_pipeline.py
LOG=/root/autodl-tmp/benchmark_1.0/filter_2_ligand_qualification_v2/runs/20260804_full_01/logs/full_screen.log

exec "$PY" "$SCRIPT" run-full --workers 24 --batch-size 200 >>"$LOG" 2>&1
