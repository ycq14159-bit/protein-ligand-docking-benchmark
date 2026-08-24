#!/usr/bin/env python3
"""High-confidence secret and large-file scan for the public repository."""

from __future__ import annotations

from pathlib import Path
import json
import re


ROOT = Path(__file__).resolve().parents[2]
PATTERNS = {
    "private_key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----"),
    "github_token": re.compile(r"\b(?:ghp|gho|ghu|ghs|ghr)_[A-Za-z0-9]{20,}\b"),
    "github_pat": re.compile(r"\bgithub_pat_[A-Za-z0-9_]{20,}\b"),
    "aws_access_key": re.compile(r"\b(?:AKIA|ASIA)[A-Z0-9]{16}\b"),
    "credential_url": re.compile(r"(?:https?|ssh)://[^\s/:]+:[^\s/@]+@"),
}


def main() -> None:
    hits = []
    large = []
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or ".git" in path.parts:
            continue
        relative = str(path.relative_to(ROOT))
        if path.stat().st_size > 20 * 1024 * 1024:
            large.append({"file": relative, "size_bytes": path.stat().st_size})
        raw = path.read_bytes()
        if b"\0" in raw[:8192]:
            continue
        text = raw.decode("utf-8", errors="replace")
        for name, pattern in PATTERNS.items():
            for match in pattern.finditer(text):
                hits.append({"file": relative, "line": text.count("\n", 0, match.start()) + 1, "pattern": name})
    result = {"validation_pass": not hits and not large, "high_confidence_secret_hits": hits, "files_over_20mb": large}
    print(json.dumps(result, indent=2))
    if not result["validation_pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
