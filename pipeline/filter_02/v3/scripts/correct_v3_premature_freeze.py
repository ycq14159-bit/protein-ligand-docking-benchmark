import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/root/autodl-tmp/benchmark_1.0/filter_2_ligand_qualification_v3")
RUN = ROOT / "runs/20260804_full_01"
AUDIT = RUN / "audit"


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


frozen = RUN / "_FROZEN.json"
current = ROOT / "CURRENT_RUN.json"
status_path = RUN / "status.json"
validation_path = RUN / "release/validation_report.json"

if not frozen.exists() or not current.exists():
    raise SystemExit("expected premature freeze markers are missing")

(AUDIT / "freeze_attempt_001_inconsistent_FROZEN.json").write_bytes(frozen.read_bytes())
(AUDIT / "freeze_attempt_001_inconsistent_CURRENT_RUN.json").write_bytes(current.read_bytes())
(AUDIT / "freeze_attempt_001_inconsistent_status.json").write_bytes(status_path.read_bytes())

status = json.loads(status_path.read_text())
validation = json.loads(validation_path.read_text())
if validation.get("validation_pass") is not True:
    raise SystemExit("release validation is not passing; refusing metadata-only correction")
if status.get("validation", {}).get("validation_pass") is not False:
    raise SystemExit("expected stale failed validation object was not found")

audit = {
    "correction_id": "v3_freeze_metadata_correction_001",
    "timestamp": datetime.now(timezone.utc).isoformat(),
    "scope": "freeze_metadata_only",
    "cause": "Successful validation retry did not replace the stale failed validation object in status.json.",
    "formal_output_rows_changed": False,
    "formal_output_files_changed": False,
    "route_semantics_changed": False,
    "placement_semantics_changed": False,
    "release_validation_sha256": sha256(validation_path),
    "action": "Archive premature markers, return run to VALIDATION_FAILED, then rerun finalize/validation with corrected status writer.",
}
(AUDIT / "v3_freeze_metadata_correction_001.json").write_text(json.dumps(audit, indent=2) + "\n")

frozen.unlink()
current.unlink()
link = ROOT / "current"
if link.is_symlink():
    link.unlink()
status["status"] = "VALIDATION_FAILED"
status.pop("frozen_at", None)
status_path.write_text(json.dumps(status, indent=2) + "\n")
print(json.dumps(audit, indent=2))
