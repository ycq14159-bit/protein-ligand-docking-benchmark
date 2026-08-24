#!/usr/bin/env python3
"""Download audited subsets of official PDB archive sources.

This script is intentionally conservative. It can build a lightweight HTTP
directory inventory, run dry-runs, and download limited smoke subsets. Full
archive runs should be launched only after reviewing audit outputs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import shutil
import sys
import tarfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from html.parser import HTMLParser
from pathlib import Path
from typing import Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin
from urllib.request import Request, urlopen

try:
    import yaml
except Exception:  # pragma: no cover
    yaml = None


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "configs" / "pdb_archive_sources_v2.yaml"
MANIFEST_FIELDS = [
    "source",
    "remote_url",
    "remote_relative_path",
    "local_path",
    "pdb_id",
    "subdirectory",
    "expected_size",
    "downloaded_size",
    "http_status",
    "resumed",
    "gzip_ok",
    "parse_ok",
    "file_type",
    "status",
    "error",
]


class LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "a":
            return
        for key, value in attrs:
            if key.lower() == "href" and value:
                self.links.append(value)


@dataclass(frozen=True)
class RemoteItem:
    source: str
    url: str
    relpath: str
    local_relpath: str
    pdb_id: str
    subdir: str
    file_type: str
    expected_size: int | None = None


def load_config(path: Path) -> dict:
    if yaml is None:
        raise SystemExit("pyyaml is required for this script")
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def fetch_bytes(url: str, timeout: int, headers: dict[str, str] | None = None) -> tuple[int, bytes, dict[str, str], str]:
    req = Request(url, headers=headers or {"User-Agent": "vs-benchmark-archive-audit/1.0"})
    try:
        with urlopen(req, timeout=timeout) as resp:
            return int(resp.status), resp.read(), dict(resp.headers), resp.geturl()
    except HTTPError as exc:
        return int(exc.code), exc.read(), dict(exc.headers), exc.geturl()


def list_http_dir(url: str, timeout: int) -> list[str]:
    status, body, headers, final_url = fetch_bytes(url, timeout)
    ctype = headers.get("Content-Type", "")
    if status >= 400:
        raise RuntimeError(f"HTTP {status} for {url}")
    if "text/html" not in ctype and b"<html" not in body[:500].lower():
        raise RuntimeError(f"Directory listing is not HTML for {url}: {ctype}")
    parser = LinkParser()
    parser.feed(body.decode("utf-8", errors="ignore"))
    links = []
    for link in parser.links:
        if link in {"../", "/"} or link.startswith("?"):
            continue
        links.append(link)
    return links


def head_size(url: str, timeout: int) -> tuple[int | None, int | None, str | None]:
    req = Request(url, method="HEAD", headers={"User-Agent": "vs-benchmark-archive-audit/1.0"})
    try:
        with urlopen(req, timeout=timeout) as resp:
            size = resp.headers.get("Content-Length")
            return int(resp.status), int(size) if size and size.isdigit() else None, resp.geturl()
    except HTTPError as exc:
        return int(exc.code), None, exc.geturl()
    except URLError:
        return None, None, None


def parse_pdb_id_from_name(source: str, name: str) -> str:
    base = Path(name).name
    if source == "mmcif":
        return base.removesuffix(".cif.gz").lower()
    if source == "pdb":
        return base.removeprefix("pdb").removesuffix(".ent.gz").lower()
    match = re.search(r"([0-9][A-Za-z0-9]{3})-pdb-bundle\.tar\.gz$", base)
    return match.group(1).lower() if match else ""


def build_inventory_for_source(source: str, cfg: dict, timeout: int, max_files: int | None) -> list[RemoteItem]:
    src = cfg["sources"][source]
    base_url = src["base_url"]
    suffix = src["expected_suffix"]
    items: list[RemoteItem] = []
    subdirs = [x.rstrip("/") for x in list_http_dir(base_url, timeout) if x.endswith("/") and len(x.rstrip("/")) == 2]
    subdirs = sorted(set(subdirs))

    if source in {"mmcif", "pdb"}:
        wanted = max_files or sys.maxsize
        per_round = []
        for subdir in subdirs:
            if len(items) >= wanted:
                break
            try:
                links = sorted(x for x in list_http_dir(urljoin(base_url, subdir + "/"), timeout) if x.endswith(suffix))
            except Exception:
                continue
            per_round.append((subdir, links))
        # Round-robin to avoid taking all files from one subdir during smoke.
        index = 0
        while len(items) < wanted:
            progressed = False
            for subdir, links in per_round:
                if index >= len(links) or len(items) >= wanted:
                    continue
                link = links[index]
                pdb_id = parse_pdb_id_from_name(source, link)
                relpath = f"{subdir}/{link}"
                url = urljoin(base_url, relpath)
                _, size, _ = head_size(url, timeout)
                local_root = "mmCIF" if source == "mmcif" else source
                local_rel = f"{local_root}/{relpath}"
                items.append(RemoteItem(source, url, relpath, local_rel, pdb_id, subdir, source, size))
                progressed = True
            if not progressed:
                break
            index += 1
        return items

    if source == "pdb_bundle":
        wanted = max_files or sys.maxsize
        for subdir in subdirs:
            if len(items) >= wanted:
                break
            try:
                entry_dirs = sorted(x.rstrip("/") for x in list_http_dir(urljoin(base_url, subdir + "/"), timeout) if x.endswith("/"))
            except Exception:
                continue
            for entry in entry_dirs:
                if len(items) >= wanted:
                    break
                try:
                    files = sorted(x for x in list_http_dir(urljoin(base_url, f"{subdir}/{entry}/"), timeout) if x.endswith(suffix))
                except Exception:
                    continue
                for fname in files[:1]:
                    pdb_id = parse_pdb_id_from_name(source, fname) or entry.lower()
                    relpath = f"{subdir}/{entry}/{fname}"
                    url = urljoin(base_url, relpath)
                    _, size, _ = head_size(url, timeout)
                    local_rel = f"pdb_bundle/{relpath}"
                    items.append(RemoteItem(source, url, relpath, local_rel, pdb_id, subdir, source, size))
                    break
        return items
    raise ValueError(source)


def ensure_free_space(output_dir: Path, min_free_gb: int) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    usage = shutil.disk_usage(output_dir)
    free_gb = usage.free / (1024**3)
    if free_gb < min_free_gb:
        raise SystemExit(f"Refusing download: free space {free_gb:.1f} GB < {min_free_gb} GB")


def verify_gzip_file(path: Path) -> tuple[bool, str]:
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1024 * 1024):
                pass
        return True, ""
    except Exception as exc:
        return False, f"{type(exc).__name__}: {exc}"


def verify_parse(item: RemoteItem, path: Path) -> tuple[bool, str]:
    try:
        if item.source == "mmcif":
            import gemmi
            with gzip.open(path, "rb") as fh:
                text = fh.read().decode("utf-8", errors="replace")
            gemmi.cif.read_string(text)
            return True, ""
        if item.source == "pdb":
            with gzip.open(path, "rt", errors="ignore") as fh:
                text = fh.read(200000)
            ok = any(line.startswith(("HEADER", "ATOM", "HETATM")) for line in text.splitlines())
            return ok, "" if ok else "No HEADER/ATOM/HETATM found in first 200 KB"
        if item.source == "pdb_bundle":
            with tarfile.open(path, "r:gz") as tf:
                names = tf.getnames()
            has_pdb = any(name.endswith(".pdb") for name in names)
            has_chain = any("chain-id" in name.lower() or "chain" in name.lower() for name in names)
            return bool(names and has_pdb), "" if has_pdb else "Bundle tar.gz has no PDB member"
    except Exception as exc:
        return False, f"{type(exc).__name__}: {exc}"
    return False, "Unknown source"


def download_one(item: RemoteItem, output_dir: Path, args: argparse.Namespace) -> dict[str, str | int | bool | None]:
    local = output_dir / item.local_relpath
    local.parent.mkdir(parents=True, exist_ok=True)
    part = local.with_suffix(local.suffix + ".part")
    row: dict[str, str | int | bool | None] = {
        "source": item.source,
        "remote_url": item.url,
        "remote_relative_path": item.relpath,
        "local_path": str(local),
        "pdb_id": item.pdb_id,
        "subdirectory": item.subdir,
        "expected_size": item.expected_size or "",
        "downloaded_size": "",
        "http_status": "",
        "resumed": False,
        "gzip_ok": "",
        "parse_ok": "",
        "file_type": item.file_type,
        "status": "",
        "error": "",
    }
    if local.exists() and local.stat().st_size > 0:
        if item.expected_size and local.stat().st_size == item.expected_size:
            row.update({"downloaded_size": local.stat().st_size, "status": "exists"})
            gzip_ok, gzerr = verify_gzip_file(local) if args.verify_gzip and local.name.endswith(".gz") else (True, "")
            parse_ok, perr = verify_parse(item, local)
            row.update({"gzip_ok": gzip_ok, "parse_ok": parse_ok, "error": gzerr or perr})
            return row
    if args.dry_run or args.inventory_only:
        row["status"] = "planned"
        return row
    headers = {"User-Agent": "vs-benchmark-archive-downloader/2.0"}
    mode = "wb"
    if args.resume and part.exists() and part.stat().st_size > 0:
        headers["Range"] = f"bytes={part.stat().st_size}-"
        mode = "ab"
        row["resumed"] = True
    last_error = ""
    for attempt in range(1, args.retries + 1):
        try:
            req = Request(item.url, headers=headers)
            with urlopen(req, timeout=args.timeout) as resp:
                row["http_status"] = int(resp.status)
                if resp.status not in {200, 206}:
                    raise RuntimeError(f"HTTP {resp.status}")
                ctype = resp.headers.get("Content-Type", "")
                if "text/html" in ctype.lower():
                    raise RuntimeError(f"Refusing to save HTML response: {ctype}")
                with part.open(mode + ("b" if "b" not in mode else "")) as out:
                    while True:
                        chunk = resp.read(1024 * 1024)
                        if not chunk:
                            break
                        out.write(chunk)
            if item.expected_size and part.stat().st_size != item.expected_size:
                raise RuntimeError(f"size mismatch: {part.stat().st_size} != {item.expected_size}")
            os.replace(part, local)
            row["downloaded_size"] = local.stat().st_size
            gzip_ok, gzerr = verify_gzip_file(local) if args.verify_gzip and local.name.endswith(".gz") else (True, "")
            parse_ok, perr = verify_parse(item, local)
            row.update({"gzip_ok": gzip_ok, "parse_ok": parse_ok})
            if gzip_ok and parse_ok:
                row["status"] = "downloaded"
            else:
                row["status"] = "qc_failed"
                row["error"] = gzerr or perr
            return row
        except Exception as exc:
            last_error = f"attempt {attempt}: {type(exc).__name__}: {exc}"
            time.sleep(min(2 * attempt, 10))
    row["status"] = "failed"
    row["error"] = last_error
    row["downloaded_size"] = part.stat().st_size if part.exists() else ""
    return row


def write_manifest(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=MANIFEST_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--source", choices=["mmcif", "pdb", "pdb_bundle", "all"], default="all")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--inventory-only", action="store_true")
    parser.add_argument("--max-files", type=int)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--verify-gzip", action="store_true")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--log-file")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg = load_config(Path(args.config))
    output_dir = Path(args.output_dir)
    min_free = int(cfg.get("download_policy", {}).get("min_free_space_gb", 200))
    ensure_free_space(output_dir, min_free)
    sources = ["mmcif", "pdb", "pdb_bundle"] if args.source == "all" else [args.source]
    all_items: list[RemoteItem] = []
    per_source_limit = args.max_files
    for source in sources:
        items = build_inventory_for_source(source, cfg, args.timeout, per_source_limit)
        all_items.extend(items)
    rows: list[dict] = []
    if args.dry_run or args.inventory_only:
        for item in all_items:
            rows.append(download_one(item, output_dir, args))
    else:
        with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
            futures = {pool.submit(download_one, item, output_dir, args): item for item in all_items}
            for fut in as_completed(futures):
                rows.append(fut.result())
                if len(rows) % 10 == 0:
                    write_manifest(Path(args.manifest), rows)
    rows.sort(key=lambda r: (str(r["source"]), str(r["remote_relative_path"])))
    write_manifest(Path(args.manifest), rows)
    summary = {
        "sources": sources,
        "planned_items": len(all_items),
        "rows": len(rows),
        "status_counts": {status: sum(1 for r in rows if r["status"] == status) for status in sorted({r["status"] for r in rows})},
        "gzip_ok_count": sum(1 for r in rows if r["gzip_ok"] is True),
        "parse_ok_count": sum(1 for r in rows if r["parse_ok"] is True),
        "manifest": str(Path(args.manifest)),
    }
    text = json.dumps(summary, indent=2, ensure_ascii=False)
    print(text)
    if args.log_file:
        Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)
        Path(args.log_file).write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
