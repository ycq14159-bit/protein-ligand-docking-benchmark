from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from audit_common import AUDIT, EXT, OUT, QC, sha256, write_tsv


def target_files() -> list[Path]:
    files = list(OUT.glob("*.csv")) + list(OUT.glob("*.tsv")) + list(QC.glob("*.tsv")) + list(QC.glob("*.json"))
    files.append(EXT / "cath" / "cath_domain_to_h.parquet")
    return sorted(p for p in files if p.name not in {"hashes.sha256", "reproducibility_qc.tsv"})


def digest_map() -> dict[str, str]:
    return {p.relative_to(AUDIT).as_posix(): sha256(p) for p in target_files() if p.exists()}


def main() -> None:
    run_all = Path(__file__).with_name("run_all.py")
    subprocess.run([sys.executable, str(run_all)], check=True)
    first = digest_map()
    subprocess.run([sys.executable, str(run_all)], check=True)
    second = digest_map()
    names = sorted(set(first) | set(second))
    rows = [{"path": name, "run1_sha256": first.get(name, "MISSING"), "run2_sha256": second.get(name, "MISSING"),
             "status": "PASS" if first.get(name) == second.get(name) else "FAIL"} for name in names]
    write_tsv(QC / "reproducibility_qc.tsv", rows)
    if not rows or any(row["status"] != "PASS" for row in rows):
        raise RuntimeError("Double-run reproducibility check failed")

    hash_targets = []
    for folder in (AUDIT / "references", AUDIT / "scripts", OUT, QC):
        hash_targets.extend(p for p in folder.rglob("*") if p.is_file() and p.name != "hashes.sha256" and "__pycache__" not in p.parts)
    lines = [f"{sha256(p)}  {p.relative_to(AUDIT).as_posix()}" for p in sorted(hash_targets)]
    (QC / "hashes.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

