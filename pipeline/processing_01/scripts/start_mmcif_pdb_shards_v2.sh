#!/usr/bin/env bash
set -euo pipefail
OUT=/root/autodl-tmp/pdb_archive_v2
INV=$OUT/inventories
LOGS=$OUT/logs
COMMON=(--recursive --no-parent --continue --timestamping --no-host-directories --reject 'index.html*' --tries=20 --timeout=60 --waitretry=5 --read-timeout=60 --retry-connrefused --progress=dot:giga)
for shard in 0 1 2; do
  wget "${COMMON[@]}" --accept '*.cif.gz' --cut-dirs=7 -P "$OUT/mmCIF" -i "$INV/mmCIF_shard${shard}_urls.txt" -o "$LOGS/mmCIF_shard${shard}.log" & echo $! > "$LOGS/mmCIF_shard${shard}.pid"
  wget "${COMMON[@]}" --accept '*.ent.gz' --cut-dirs=7 -P "$OUT/pdb" -i "$INV/pdb_shard${shard}_urls.txt" -o "$LOGS/pdb_shard${shard}.log" & echo $! > "$LOGS/pdb_shard${shard}.pid"
done
wait
