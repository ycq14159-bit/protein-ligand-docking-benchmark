#!/usr/bin/env bash
set -euo pipefail
source /root/miniconda3/etc/profile.d/conda.sh
ROOT=/root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality_v2
RUN=$ROOT/runs/20260814_full_01
CONFIG=$RUN/config_snapshot.yaml
conda activate entry-work-packages-v1
python "$ROOT/scripts/filter3_v2_pipeline.py" --config "$CONFIG" --run-dir "$RUN" --workers 8
conda activate posebusters-audit
PYTHONPATH=/root/miniconda3/envs/entry-work-packages-v1/lib/python3.10/site-packages \
python "$ROOT/scripts/filter3_v2_posebusters.py" --config "$CONFIG" --run-dir "$RUN" --workers 8
conda activate entry-work-packages-v1
python "$ROOT/scripts/filter3_v2_apply_posebusters.py" --config "$CONFIG" --run-dir "$RUN" --workers 8
python "$ROOT/scripts/filter3_v2_finalize.py" --run-dir "$RUN"
python "$ROOT/scripts/validate_filter3_v2_release.py" --run-dir "$RUN"
