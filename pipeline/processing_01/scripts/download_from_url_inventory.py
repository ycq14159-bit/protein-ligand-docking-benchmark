#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import os
import queue
import shutil
import sys
import threading
import time
from datetime import datetime
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

FIELDS = [
    "timestamp",
    "source",
    "url",
    "pdb_id",
    "local_path",
    "expected_size",
    "downloaded_size",
    "resumed",
    "attempts",
    "http_status",
    "gzip_ok",
    "status",
    "error",
]


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def gzip_ok(path: Path) -> bool:
    try:
        if not path.exists() or path.stat().st_size == 0:
            return False
        with gzip.open(path, "rb") as fh:
            while fh.read(1024 * 1024):
                pass
        return True
    except Exception:
        return False


def looks_html(path: Path) -> bool:
    try:
        with path.open("rb") as fh:
            head = fh.read(256).lstrip().lower()
        return head.startswith(b"<!doctype html") or head.startswith(b"<html")
    except Exception:
        return False


def parse_expected_size(status: int, headers: dict[str, str], existing_size: int) -> int | str:
    if status == 206:
        cr = headers.get("Content-Range", "")
        if "/" in cr:
            total = cr.rsplit("/", 1)[-1]
            return int(total) if total.isdigit() else ""
    cl = headers.get("Content-Length", "")
    if cl.isdigit():
        return int(cl) + existing_size if status == 206 else int(cl)
    return ""


def row_base(plan_row: dict[str, str]) -> dict[str, object]:
    return {
        "timestamp": now(),
        "source": plan_row["source"],
        "url": plan_row["url"],
        "pdb_id": plan_row.get("pdb_id", ""),
        "local_path": plan_row["local_path"],
        "expected_size": "",
        "downloaded_size": "",
        "resumed": "false",
        "attempts": 0,
        "http_status": "",
        "gzip_ok": "",
        "status": "",
        "error": "",
    }


def disk_free_gb(path: Path) -> float:
    path.mkdir(parents=True, exist_ok=True)
    return shutil.disk_usage(path).free / (1024 ** 3)


class Downloader:
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.output_root = Path(args.output_root)
        self.stop_event = threading.Event()

    def should_stop_for_disk(self) -> bool:
        return disk_free_gb(self.output_root) < self.args.min_free_gb

    def download_one(self, plan_row: dict[str, str]) -> dict[str, object]:
        out = row_base(plan_row)
        local = Path(plan_row["local_path"])
        part = Path(str(local) + ".part")
        local.parent.mkdir(parents=True, exist_ok=True)

        if local.exists():
            ok = gzip_ok(local) if self.args.verify_gzip else local.stat().st_size > 0
            out.update({
                "downloaded_size": local.stat().st_size,
                "gzip_ok": str(ok).lower(),
                "status": "skipped_existing" if ok else "corrupt_existing",
                "error": "" if ok else "existing final file failed gzip/nonzero check",
            })
            return out

        if plan_row.get("action") == "review_corrupt_existing":
            out.update({"status": "corrupt_existing", "error": "marked corrupt in download plan"})
            return out

        if self.should_stop_for_disk():
            self.stop_event.set()
            out.update({"status": "disk_guard_stop", "error": f"free space below {self.args.min_free_gb} GB"})
            return out

        last_error = ""
        for attempt in range(1, self.args.retries + 1):
            if self.stop_event.is_set():
                out.update({"attempts": attempt - 1, "status": "disk_guard_stop", "error": "stop event set"})
                return out
            existing = part.stat().st_size if part.exists() else 0
            headers = {"User-Agent": "vs-benchmark-url-inventory-downloader/3.0"}
            mode = "wb"
            resumed = False
            if self.args.resume and existing > 0:
                headers["Range"] = f"bytes={existing}-"
                mode = "ab"
                resumed = True
            try:
                req = Request(plan_row["url"], headers=headers)
                with urlopen(req, timeout=self.args.timeout) as resp:
                    status = int(resp.status)
                    out["http_status"] = status
                    ctype = resp.headers.get("Content-Type", "")
                    if "text/html" in ctype.lower():
                        raise RuntimeError(f"refusing HTML content-type: {ctype}")
                    if resumed and status != 206:
                        raise RuntimeError(f"server did not honor Range for existing part: HTTP {status}")
                    if status not in (200, 206):
                        raise RuntimeError(f"HTTP {status}")
                    out["expected_size"] = parse_expected_size(status, dict(resp.headers), existing)
                    with part.open(mode) as fh:
                        while True:
                            chunk = resp.read(1024 * 1024)
                            if not chunk:
                                break
                            fh.write(chunk)
                if not part.exists() or part.stat().st_size == 0:
                    raise RuntimeError("downloaded part is empty")
                if looks_html(part):
                    raise RuntimeError("downloaded part looks like HTML")
                ok = gzip_ok(part) if self.args.verify_gzip else True
                if not ok:
                    raise RuntimeError("gzip verification failed")
                if local.exists():
                    raise RuntimeError("final path appeared before atomic rename")
                os.replace(part, local)
                out.update({
                    "downloaded_size": local.stat().st_size,
                    "resumed": str(resumed).lower(),
                    "attempts": attempt,
                    "gzip_ok": "true",
                    "status": "resumed_completed" if resumed else "downloaded",
                    "error": "",
                })
                return out
            except Exception as exc:
                last_error = f"{type(exc).__name__}: {exc}"
                time.sleep(min(2 * attempt, 10))
        out.update({
            "downloaded_size": part.stat().st_size if part.exists() else "",
            "resumed": str(part.exists()).lower(),
            "attempts": self.args.retries,
            "gzip_ok": "false",
            "status": "failed",
            "error": last_error,
        })
        return out


def read_plan(path: Path, source: str) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = [r for r in reader if not source or r.get("source") == source]
    seen: set[str] = set()
    deduped: list[dict[str, str]] = []
    duplicates = 0
    for row in rows:
        lp = row["local_path"]
        if lp in seen:
            duplicates += 1
            continue
        seen.add(lp)
        deduped.append(row)
    if duplicates:
        raise SystemExit(f"Refusing to run: duplicate local_path count={duplicates}")
    return deduped


def write_row(writer: csv.DictWriter, fh, row: dict[str, object], count: int) -> None:
    writer.writerow(row)
    if count % 100 == 0:
        fh.flush()
        os.fsync(fh.fileno())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True)
    ap.add_argument("--inventory", required=True)
    ap.add_argument("--output-root", required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument("--timeout", type=int, default=60)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--log-file", required=True)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--verify-gzip", action="store_true")
    ap.add_argument("--min-free-gb", type=int, default=200)
    args = ap.parse_args()

    rows = read_plan(Path(args.inventory), args.source)
    Path(args.manifest).parent.mkdir(parents=True, exist_ok=True)
    Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)
    dl = Downloader(args)

    pending = queue.Queue()
    immediate: list[dict[str, object]] = []
    for row in rows:
        if row.get("action") == "skip_completed":
            out = row_base(row)
            local = Path(row["local_path"])
            out.update({
                "downloaded_size": local.stat().st_size if local.exists() else "",
                "gzip_ok": "true",
                "status": "skipped_existing",
            })
            immediate.append(out)
        elif row.get("action") == "review_corrupt_existing":
            out = row_base(row)
            out.update({"status": "corrupt_existing", "error": "marked corrupt in download plan"})
            immediate.append(out)
        else:
            pending.put(row)

    result_q: queue.Queue[dict[str, object]] = queue.Queue()

    def worker() -> None:
        while not dl.stop_event.is_set():
            try:
                item = pending.get_nowait()
            except queue.Empty:
                return
            try:
                result_q.put(dl.download_one(item))
            finally:
                pending.task_done()

    threads = [threading.Thread(target=worker, daemon=True) for _ in range(max(1, args.workers))]
    start = time.time()
    with Path(args.manifest).open("w", newline="", encoding="utf-8") as mf, Path(args.log_file).open("w", encoding="utf-8") as lf:
        writer = csv.DictWriter(mf, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        count = 0
        status_counts: dict[str, int] = {}
        for row in immediate:
            count += 1
            status_counts[str(row["status"])] = status_counts.get(str(row["status"]), 0) + 1
            write_row(writer, mf, row, count)
        for t in threads:
            t.start()
        while any(t.is_alive() for t in threads) or not result_q.empty():
            try:
                row = result_q.get(timeout=2)
            except queue.Empty:
                continue
            count += 1
            status_counts[str(row["status"])] = status_counts.get(str(row["status"]), 0) + 1
            write_row(writer, mf, row, count)
            if count % 100 == 0:
                elapsed = max(time.time() - start, 1)
                msg = f"{now()} source={args.source} processed={count}/{len(rows)} pending={pending.qsize()} rate={count/elapsed:.2f}/s status={status_counts}"
                print(msg, flush=True)
                lf.write(msg + "\n")
                lf.flush()
        for t in threads:
            t.join()
        final = f"{now()} source={args.source} DONE processed={count}/{len(rows)} status={status_counts}"
        print(final, flush=True)
        lf.write(final + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
