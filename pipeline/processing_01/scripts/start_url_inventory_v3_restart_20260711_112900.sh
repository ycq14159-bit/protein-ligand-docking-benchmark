#!/usr/bin/env bash
set -euo pipefail
ROOT=/root/autodl-tmp/pdb_archive_v2
SCRIPT=$ROOT/scripts/download_from_url_inventory.py
LOGDIR=/root/autodl-tmp/pdb_archive_v2/logs/url_inventory_v3_restart_20260711_112900
mkdir -p "$LOGDIR"
CMD1="/root/miniconda3/bin/python $SCRIPT --source mmcif --inventory $ROOT/manifests/mmcif_download_plan.tsv --output-root $ROOT --workers 4 --retries 3 --timeout 60 --manifest $ROOT/manifests/mmcif_download_runtime_restart_20260711_112900.tsv --log-file $LOGDIR/mmcif.log --resume --verify-gzip"
CMD2="/root/miniconda3/bin/python $SCRIPT --source pdb --inventory $ROOT/manifests/pdb_download_plan.tsv --output-root $ROOT --workers 4 --retries 3 --timeout 60 --manifest $ROOT/manifests/pdb_download_runtime_restart_20260711_112900.tsv --log-file $LOGDIR/pdb.log --resume --verify-gzip"
{
  echo "started_at=$(date '+%F %T %Z')"
  echo "session=pdb_archive_url_inventory_v3_restart"
  echo "mmcif: $CMD1"
  echo "pdb: $CMD2"
} > "$LOGDIR/start_commands.txt"
set +e
bash -lc "$CMD1" > "$LOGDIR/mmcif.stdout.log" 2> "$LOGDIR/mmcif.stderr.log" & echo $! > "$LOGDIR/mmcif.pid"
bash -lc "$CMD2" > "$LOGDIR/pdb.stdout.log" 2> "$LOGDIR/pdb.stderr.log" & echo $! > "$LOGDIR/pdb.pid"
wait
code=$?
echo "download_children_exited_at=$(date '+%F %T %Z') exit_code=$code" | tee -a "$LOGDIR/supervisor.log"
echo "Keeping screen session alive for diagnosis. Press Ctrl-C or exit to close." | tee -a "$LOGDIR/supervisor.log"
while true; do sleep 3600; done
