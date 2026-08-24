#!/usr/bin/env bash
set -euo pipefail
OUT=/root/autodl-tmp/pdb_archive_v2
mkdir -p "$OUT"/{mmCIF,pdb,pdb_bundle,inventories,manifests,logs,failed,scripts}
COMMON=(--recursive --no-parent --continue --timestamping --no-host-directories --reject 'index.html*' --tries=20 --timeout=60 --waitretry=5 --read-timeout=60 --retry-connrefused --progress=dot:giga)
# Keep each source separate. Wget recursively traverses all official subdirectories and skips complete files via -c/-N.
wget "${COMMON[@]}" --accept '*.cif.gz' --cut-dirs=7 -P "$OUT/mmCIF" https://files.rcsb.org/pub/pdb/data/structures/divided/mmCIF/ -o "$OUT/logs/mmCIF_wget.log" & echo $! > "$OUT/logs/mmCIF.pid"
wget "${COMMON[@]}" --accept '*.ent.gz' --cut-dirs=7 -P "$OUT/pdb" https://files.rcsb.org/pub/pdb/data/structures/divided/pdb/ -o "$OUT/logs/pdb_wget.log" & echo $! > "$OUT/logs/pdb.pid"
wget "${COMMON[@]}" --accept '*-pdb-bundle.tar.gz' --cut-dirs=5 -P "$OUT/pdb_bundle" https://files.wwpdb.org/pub/pdb/compatible/pdb_bundle/ -o "$OUT/logs/pdb_bundle_wget.log" & echo $! > "$OUT/logs/pdb_bundle.pid"
wait
