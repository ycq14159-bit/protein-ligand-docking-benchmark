#!/usr/bin/env bash
set -u
ST=/root/autodl-tmp/benchmark_1.0/processing_03_direct_contact_qualification
AUDIT="$ST/plip_validation"
RUN="$AUDIT/runs/20260812_full_01"
PY=/root/miniconda3/envs/entry-work-packages-v1/bin/python
LOG="$RUN/logs/health_monitor.log"

while true; do
  now=$(date -Is)
  status=$("$PY" -c "import json; print(json.load(open('$RUN/status.json')).get('status','MISSING'))" 2>/dev/null || echo MISSING)
  pairs=$("$PY" -c "import json; print(json.load(open('$RUN/status.json')).get('pairs_this_attempt',0))" 2>/dev/null || echo 0)
  checkpoints=$(find "$RUN/work/checkpoints" -type f 2>/dev/null | wc -l)
  plip_procs=$(pgrep -fc '/root/miniconda3/envs/plip-audit/bin/plip' 2>/dev/null || true)
  mem=$(awk '/MemAvailable:/{print $2}' /proc/meminfo)
  disk=$(df -Pk /root/autodl-tmp | awk 'NR==2{print $4}')
  echo "$now status=$status pairs=$pairs checkpoints=$checkpoints plip_procs=$plip_procs mem_available_kb=$mem disk_free_kb=$disk" >> "$LOG"
  if [[ "$status" == "COMPLETED" ]]; then
    echo "$now COMPLETED starting_finalize" >> "$LOG"
    "$PY" "$AUDIT/scripts/finalize_processing3_plip.py" >> "$RUN/logs/finalize.log" 2>&1
    rc=$?
    echo "$(date -Is) FINALIZE_EXIT=$rc" >> "$LOG"
    exit "$rc"
  fi
  if ! screen -ls 2>/dev/null | grep -q '[.]processing3_plip_full'; then
    echo "$now ALERT runner_screen_missing status=$status" >> "$LOG"
  fi
  sleep 300
done
