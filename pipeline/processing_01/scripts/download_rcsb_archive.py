from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


PROJECT_ROOT = Path(__file__).resolve().parents[1]

MMCIF_BASE_URL = "https://files.rcsb.org/pub/pdb/data/structures/divided/mmCIF"
PDB_BASE_URL = "https://files.rcsb.org/pub/pdb/data/structures/divided/pdb"
BUNDLE_BASE_URL = "https://files.wwpdb.org/pub/pdb/compatible/pdb_bundle"

REQUIRED_CANDIDATE_COLUMNS = [
    "pdb_id",
    "resolution",
    "release_date",
    "experimental_method",
    "organism",
    "candidate_notes",
]
MANIFEST_COLUMNS = [
    "pdb_id",
    "preferred_format",
    "downloaded_format",
    "url",
    "local_path",
    "status",
    "error",
    "file_size_bytes",
]


@dataclass
class DownloadResult:
    pdb_id: str
    preferred_format: str
    downloaded_format: str
    url: str
    local_path: Path | None
    status: str
    error: str
    file_size_bytes: int


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download selected RCSB structures from official divided PDB archives."
    )
    parser.add_argument(
        "--candidate_csv",
        default="data/rcsb_archive_raw/manifests/rcsb_candidate_ids.csv",
        help="CSV with candidate pdb_id values and metadata.",
    )
    parser.add_argument(
        "--out_dir",
        default="data/rcsb_archive_raw",
        help="Root output directory for downloaded structures and manifests.",
    )
    parser.add_argument(
        "--prefer",
        choices=["mmcif", "pdb", "pdb_bundle"],
        default="mmcif",
        help="Preferred download format.",
    )
    parser.add_argument(
        "--fallback",
        default="pdb,pdb_bundle",
        help="Comma-separated fallback formats, for example: pdb or pdb,pdb_bundle.",
    )
    parser.add_argument(
        "--max_samples",
        type=int,
        default=0,
        help="Maximum candidate rows to try; use 0 for all rows.",
    )
    parser.add_argument("--timeout", type=int, default=60, help="HTTP timeout in seconds.")
    parser.add_argument("--retries", type=int, default=2, help="Retries per URL after the first attempt.")
    parser.add_argument("--progress_every", type=int, default=100, help="Print progress every N candidates.")
    args = parser.parse_args()

    candidate_csv = resolve_project_path(args.candidate_csv)
    out_dir = resolve_project_path(args.out_dir)
    manifests_dir = out_dir / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)

    if not candidate_csv.exists():
        example_csv = manifests_dir / "rcsb_candidate_ids.example.csv"
        write_example_candidate_csv(example_csv)
        if candidate_csv.name == "rcsb_candidate_ids.csv":
            shutil.copyfile(example_csv, candidate_csv)
            print(f"Candidate CSV was missing; wrote example candidates to: {candidate_csv}")
        else:
            print(f"Candidate CSV was missing; wrote example file to: {example_csv}")
            print("Re-run with --candidate_csv pointing to an existing candidate CSV.")
            return 2

    all_candidates = read_candidates(candidate_csv)
    candidates = list(all_candidates)
    if args.max_samples and args.max_samples > 0:
        candidates = candidates[: args.max_samples]
    if not candidates:
        raise ValueError(f"No candidate rows found in {candidate_csv}")

    formats = build_format_order(args.prefer, args.fallback)
    ensure_output_dirs(out_dir)

    results = []
    for index, row in enumerate(candidates, start=1):
        result = download_one(row["pdb_id"], formats, args.prefer, out_dir, timeout=args.timeout, retries=args.retries)
        results.append(result)
        if args.progress_every and (index % args.progress_every == 0 or index == len(candidates)):
            print(
                f"[{index}/{len(candidates)}] downloaded={sum(1 for r in results if r.status == 'downloaded')} "
                f"exists={sum(1 for r in results if r.status == 'exists')} "
                f"failed={sum(1 for r in results if r.status == 'failed')}",
                flush=True,
            )

    manifest_path = manifests_dir / "download_manifest.csv"
    write_manifest(manifest_path, results)
    report_path = manifests_dir / "download_report.md"
    write_report(report_path, candidate_csv, out_dir, formats, results, candidate_total=len(all_candidates))

    print(f"Tried candidates: {len(results)}")
    print(f"Downloaded: {sum(1 for r in results if r.status == 'downloaded')}")
    print(f"Already existed: {sum(1 for r in results if r.status == 'exists')}")
    print(f"Failed: {sum(1 for r in results if r.status == 'failed')}")
    print(f"Wrote manifest: {manifest_path}")
    print(f"Wrote report: {report_path}")
    return 0 if all(r.status in {"downloaded", "exists"} for r in results) else 1


def resolve_project_path(value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return PROJECT_ROOT / path


def write_example_candidate_csv(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        ("1a4w", "1.80", "1998-04-15", "X-ray diffraction", "Homo sapiens", "Example smoke-test candidate"),
        ("1bcu", "2.00", "1998-08-26", "X-ray diffraction", "Escherichia coli", "Example smoke-test candidate"),
        ("1cbs", "2.00", "1995-01-26", "X-ray diffraction", "Homo sapiens", "Example smoke-test candidate"),
        ("1hvr", "2.00", "1994-01-31", "X-ray diffraction", "HIV-1", "Example smoke-test candidate"),
        ("1m17", "2.60", "2002-08-07", "X-ray diffraction", "Homo sapiens", "Example smoke-test candidate"),
        ("1qcf", "2.30", "1999-04-06", "X-ray diffraction", "Homo sapiens", "Example smoke-test candidate"),
        ("1stp", "2.60", "1992-10-15", "X-ray diffraction", "Streptomyces avidinii", "Example smoke-test candidate"),
        ("2br1", "1.90", "2005-01-25", "X-ray diffraction", "Homo sapiens", "Example smoke-test candidate"),
        ("2hyy", "1.90", "2006-08-29", "X-ray diffraction", "Homo sapiens", "Example smoke-test candidate"),
        ("2qbr", "2.20", "2007-07-17", "X-ray diffraction", "Homo sapiens", "Example smoke-test candidate"),
        ("3ert", "1.90", "1998-12-09", "X-ray diffraction", "Homo sapiens", "Example smoke-test candidate"),
        ("3pbl", "2.00", "2010-11-24", "X-ray diffraction", "Homo sapiens", "Example smoke-test candidate"),
        ("4hjo", "2.40", "2013-02-13", "X-ray diffraction", "Homo sapiens", "Example smoke-test candidate"),
        ("4j8m", "1.80", "2013-03-27", "X-ray diffraction", "Homo sapiens", "Example smoke-test candidate"),
        ("4tmn", "1.70", "1994-07-31", "X-ray diffraction", "Bacillus thermoproteolyticus", "Example smoke-test candidate"),
        ("5hvp", "2.00", "1996-05-31", "X-ray diffraction", "HIV-1", "Example smoke-test candidate"),
        ("5l2s", "1.85", "2016-09-14", "X-ray diffraction", "Homo sapiens", "Example smoke-test candidate"),
        ("6lu7", "2.16", "2020-02-05", "X-ray diffraction", "SARS-CoV-2", "Example smoke-test candidate"),
        ("7jrn", "1.75", "2020-11-11", "X-ray diffraction", "SARS-CoV-2", "Example smoke-test candidate"),
        ("8g6j", "1.80", "2023-04-19", "X-ray diffraction", "Homo sapiens", "Example smoke-test candidate"),
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(REQUIRED_CANDIDATE_COLUMNS)
        writer.writerows(rows)


def read_candidates(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = [col for col in REQUIRED_CANDIDATE_COLUMNS if col not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"Candidate CSV is missing required columns: {', '.join(missing)}")
        rows = []
        seen = set()
        for row in reader:
            pdb_id = normalize_pdb_id(row.get("pdb_id", ""))
            if not pdb_id or pdb_id in seen:
                continue
            row["pdb_id"] = pdb_id
            rows.append(row)
            seen.add(pdb_id)
    return rows


def normalize_pdb_id(value: str) -> str:
    pdb_id = value.strip().lower()
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        return ""
    return pdb_id


def build_format_order(prefer: str, fallback: str) -> list[str]:
    formats = [prefer]
    for item in fallback.split(","):
        fmt = item.strip().lower()
        if not fmt:
            continue
        if fmt not in {"mmcif", "pdb", "pdb_bundle"}:
            raise ValueError(f"Unsupported fallback format: {fmt}")
        if fmt not in formats:
            formats.append(fmt)
    return formats


def ensure_output_dirs(out_dir: Path) -> None:
    for rel in ["mmcif", "pdb", "pdb_bundle", "manifests"]:
        (out_dir / rel).mkdir(parents=True, exist_ok=True)


def divided_subdir(pdb_id: str) -> str:
    return pdb_id[1:3]


def format_target(format_name: str, pdb_id: str, out_dir: Path) -> tuple[str, Path]:
    subdir = divided_subdir(pdb_id)
    if format_name == "mmcif":
        return (
            f"{MMCIF_BASE_URL}/{subdir}/{pdb_id}.cif.gz",
            out_dir / "mmcif" / f"{pdb_id}.cif.gz",
        )
    if format_name == "pdb":
        return (
            f"{PDB_BASE_URL}/{subdir}/pdb{pdb_id}.ent.gz",
            out_dir / "pdb" / f"pdb{pdb_id}.ent.gz",
        )
    if format_name == "pdb_bundle":
        return (
            f"{BUNDLE_BASE_URL}/{subdir}/{pdb_id}/{pdb_id}-pdb-bundle.tar.gz",
            out_dir / "pdb_bundle" / f"{pdb_id}-pdb-bundle.tar.gz",
        )
    raise ValueError(f"Unsupported format: {format_name}")


def download_one(
    pdb_id: str,
    formats: Iterable[str],
    preferred_format: str,
    out_dir: Path,
    timeout: int,
    retries: int,
) -> DownloadResult:
    errors = []
    attempted_urls = []
    for format_name in formats:
        url, local_path = format_target(format_name, pdb_id, out_dir)
        attempted_urls.append(f"{format_name}={url}")
        print(f"[{pdb_id}] trying {format_name}: {url}", flush=True)
        if local_path.exists() and local_path.stat().st_size > 0:
            print(f"[{pdb_id}] exists: {local_path}", flush=True)
            return DownloadResult(
                pdb_id=pdb_id,
                preferred_format=preferred_format,
                downloaded_format=format_name,
                url=url,
                local_path=local_path,
                status="exists",
                error="",
                file_size_bytes=local_path.stat().st_size,
            )
        try:
            fetch_with_resume(url, local_path, timeout=timeout, retries=retries)
            print(f"[{pdb_id}] downloaded {format_name}: {local_path}", flush=True)
            return DownloadResult(
                pdb_id=pdb_id,
                preferred_format=preferred_format,
                downloaded_format=format_name,
                url=url,
                local_path=local_path,
                status="downloaded",
                error="",
                file_size_bytes=local_path.stat().st_size,
            )
        except Exception as exc:
            print(f"[{pdb_id}] failed {format_name}: {exc}", flush=True)
            errors.append(f"{format_name}: {exc}")
    return DownloadResult(
        pdb_id=pdb_id,
        preferred_format=preferred_format,
        downloaded_format="",
        url="; ".join(attempted_urls),
        local_path=None,
        status="failed",
        error=" | ".join(errors),
        file_size_bytes=0,
    )


def fetch_with_resume(url: str, local_path: Path, timeout: int, retries: int) -> None:
    local_path.parent.mkdir(parents=True, exist_ok=True)
    part_path = local_path.with_suffix(local_path.suffix + ".part")
    last_error: Exception | None = None
    for attempt in range(retries + 1):
        try:
            resume_from = part_path.stat().st_size if part_path.exists() else 0
            request = Request(url, headers={"User-Agent": "vs-benchmark-rcsb-downloader/0.1"})
            if resume_from:
                request.add_header("Range", f"bytes={resume_from}-")
            with urlopen(request, timeout=timeout) as response:
                if resume_from and getattr(response, "status", None) == 200:
                    resume_from = 0
                    part_path.unlink(missing_ok=True)
                mode = "ab" if resume_from else "wb"
                with part_path.open(mode) as out_handle:
                    while True:
                        chunk = response.read(1024 * 256)
                        if not chunk:
                            break
                        out_handle.write(chunk)
            if part_path.stat().st_size <= 0:
                raise RuntimeError("downloaded file is empty")
            validate_gzip_if_needed(part_path)
            part_path.replace(local_path)
            return
        except HTTPError as exc:
            last_error = exc
            if exc.code == 404:
                break
        except (URLError, TimeoutError, OSError, RuntimeError) as exc:
            last_error = exc
        if attempt < retries:
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(str(last_error) if last_error else "download failed")


def validate_gzip_if_needed(path: Path) -> None:
    if not path.name.endswith(".gz"):
        return
    with gzip.open(path, "rb") as handle:
        handle.read(1)


def write_manifest(path: Path, results: list[DownloadResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS)
        writer.writeheader()
        for result in results:
            writer.writerow(
                {
                    "pdb_id": result.pdb_id,
                    "preferred_format": result.preferred_format,
                    "downloaded_format": result.downloaded_format,
                    "url": result.url,
                    "local_path": str(result.local_path) if result.local_path else "",
                    "status": result.status,
                    "error": result.error,
                    "file_size_bytes": result.file_size_bytes,
                }
            )


def write_report(
    path: Path,
    candidate_csv: Path,
    out_dir: Path,
    formats: list[str],
    results: list[DownloadResult],
    candidate_total: int,
) -> None:
    counts = {status: sum(1 for r in results if r.status == status) for status in ["downloaded", "exists", "failed"]}
    mmcif_dir = out_dir / "mmcif"
    local_mmcif_files = list(mmcif_dir.glob("*.cif.gz")) if mmcif_dir.exists() else []
    local_mmcif_size = sum(path.stat().st_size for path in local_mmcif_files if path.is_file())
    by_format: dict[str, int] = {}
    by_format_downloaded: dict[str, int] = {}
    by_format_exists: dict[str, int] = {}
    for result in results:
        if result.downloaded_format:
            by_format[result.downloaded_format] = by_format.get(result.downloaded_format, 0) + 1
            if result.status == "downloaded":
                by_format_downloaded[result.downloaded_format] = by_format_downloaded.get(result.downloaded_format, 0) + 1
            if result.status == "exists":
                by_format_exists[result.downloaded_format] = by_format_exists.get(result.downloaded_format, 0) + 1

    lines = [
        "# RCSB archive download report",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        f"Candidate CSV: `{candidate_csv}`",
        f"Output directory: `{out_dir}`",
        f"Format order: `{', '.join(formats)}`",
        "",
        "## Summary",
        "",
        f"- Candidate total: {candidate_total}",
        f"- Tried candidates: {len(results)}",
        f"- Downloaded: {counts['downloaded']}",
        f"- Already existed: {counts['exists']}",
        f"- Failed: {counts['failed']}",
        f"- mmCIF success: {by_format.get('mmcif', 0)}",
        f"- PDB fallback success: {by_format.get('pdb', 0)}",
        f"- pdb_bundle fallback success: {by_format.get('pdb_bundle', 0)}",
        f"- Local mmCIF file count: {len(local_mmcif_files)}",
        f"- Local mmCIF total size MB: {local_mmcif_size / (1024 * 1024):.2f}",
        "",
        "## Formats",
        "",
    ]
    if by_format:
        lines.extend(f"- {fmt}: {count}" for fmt, count in sorted(by_format.items()))
    else:
        lines.append("- None")
    lines.extend(["", "## Downloaded by format", ""])
    if by_format_downloaded:
        lines.extend(f"- {fmt}: {count}" for fmt, count in sorted(by_format_downloaded.items()))
    else:
        lines.append("- None")
    lines.extend(["", "## Already existed by format", ""])
    if by_format_exists:
        lines.extend(f"- {fmt}: {count}" for fmt, count in sorted(by_format_exists.items()))
    else:
        lines.append("- None")
    lines.extend(["", "## Failed entries", ""])
    failed = [r for r in results if r.status == "failed"]
    if failed:
        lines.extend(f"- {r.pdb_id}: {r.error}" for r in failed)
    else:
        lines.append("- None")
    lines.extend(["", "## Files", ""])
    for result in results:
        lines.append(
            f"- {result.pdb_id}: {result.status}, {result.downloaded_format or 'none'}, "
            f"{result.file_size_bytes} bytes, `{result.local_path or ''}`"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
