from __future__ import annotations

import argparse
import csv
import json
import shutil
import socket
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, urlopen


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MMCIF_BASE_URL = "https://files.rcsb.org/pub/pdb/data/structures/divided/mmCIF"
PDB_BASE_URL = "https://files.rcsb.org/pub/pdb/data/structures/divided/pdb"


@dataclass
class MethodResult:
    method: str
    success: bool
    status_code: str
    error: str


@dataclass
class UrlDiagnostics:
    pdb_id: str
    format_name: str
    url: str
    host: str
    dns_success: bool
    dns_addresses: str
    dns_error: str
    method_results: list[MethodResult]
    likely_issue: str


def main() -> int:
    parser = argparse.ArgumentParser(description="Diagnose single-file RCSB archive download connectivity.")
    parser.add_argument(
        "--candidate_csv",
        default="data/rcsb_archive_raw/manifests/rcsb_candidate_ids.csv",
        help="Candidate CSV containing pdb_id.",
    )
    parser.add_argument(
        "--report",
        default="data/rcsb_archive_raw/manifests/download_diagnostics_report.md",
        help="Markdown diagnostics report path.",
    )
    parser.add_argument("--max_ids", type=int, default=3, help="Number of candidate PDB IDs to diagnose.")
    parser.add_argument("--timeout", type=int, default=20, help="Per-method timeout in seconds.")
    args = parser.parse_args()

    candidate_csv = resolve_project_path(args.candidate_csv)
    report_path = resolve_project_path(args.report)
    pdb_ids = read_pdb_ids(candidate_csv, args.max_ids)
    if not pdb_ids:
        raise ValueError(f"No pdb_id rows found in {candidate_csv}")

    curl_path = find_curl()
    diagnostics: list[UrlDiagnostics] = []
    for pdb_id in pdb_ids:
        for format_name, url in build_urls(pdb_id):
            diagnostics.append(diagnose_url(pdb_id, format_name, url, timeout=args.timeout, curl_path=curl_path))

    write_report(report_path, candidate_csv, diagnostics, curl_path)
    print(f"Diagnosed {len(pdb_ids)} PDB IDs and {len(diagnostics)} URLs")
    print(f"Wrote report: {report_path}")
    return 0


def resolve_project_path(value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return PROJECT_ROOT / path


def read_pdb_ids(path: Path, max_ids: int) -> list[str]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if "pdb_id" not in (reader.fieldnames or []):
            raise ValueError(f"Missing pdb_id column in {path}")
        pdb_ids = []
        seen = set()
        for row in reader:
            pdb_id = row.get("pdb_id", "").strip().lower()
            if len(pdb_id) != 4 or not pdb_id.isalnum() or pdb_id in seen:
                continue
            pdb_ids.append(pdb_id)
            seen.add(pdb_id)
            if max_ids > 0 and len(pdb_ids) >= max_ids:
                break
    return pdb_ids


def build_urls(pdb_id: str) -> list[tuple[str, str]]:
    subdir = pdb_id[1:3]
    return [
        ("mmcif", f"{MMCIF_BASE_URL}/{subdir}/{pdb_id}.cif.gz"),
        ("pdb", f"{PDB_BASE_URL}/{subdir}/pdb{pdb_id}.ent.gz"),
    ]


def diagnose_url(pdb_id: str, format_name: str, url: str, timeout: int, curl_path: str | None) -> UrlDiagnostics:
    host = urlparse(url).hostname or ""
    dns_success, dns_addresses, dns_error = resolve_dns(host)
    results = [
        test_urllib(url, timeout),
        test_requests(url, timeout),
        test_powershell(url, timeout),
    ]
    if curl_path:
        results.append(test_curl(url, timeout, curl_path))
    else:
        results.append(MethodResult("curl -I", False, "", "curl executable not found"))
    likely_issue = classify_issue(dns_success, results)
    return UrlDiagnostics(
        pdb_id=pdb_id,
        format_name=format_name,
        url=url,
        host=host,
        dns_success=dns_success,
        dns_addresses=dns_addresses,
        dns_error=dns_error,
        method_results=results,
        likely_issue=likely_issue,
    )


def resolve_dns(host: str) -> tuple[bool, str, str]:
    try:
        records = socket.getaddrinfo(host, 443, type=socket.SOCK_STREAM)
        addresses = sorted({record[4][0] for record in records})
        return True, ", ".join(addresses), ""
    except OSError as exc:
        return False, "", str(exc)


def test_urllib(url: str, timeout: int) -> MethodResult:
    try:
        request = Request(url, headers={"User-Agent": "vs-benchmark-rcsb-diagnostics/0.1", "Range": "bytes=0-0"})
        with urlopen(request, timeout=timeout) as response:
            response.read(1)
            return MethodResult("Python urllib.request.urlopen", True, str(getattr(response, "status", "")), "")
    except HTTPError as exc:
        return MethodResult("Python urllib.request.urlopen", False, str(exc.code), str(exc))
    except (URLError, TimeoutError, OSError) as exc:
        return MethodResult("Python urllib.request.urlopen", False, "", str(exc))


def test_requests(url: str, timeout: int) -> MethodResult:
    try:
        import requests
    except ImportError as exc:
        return MethodResult("Python requests.get", False, "", f"requests unavailable: {exc}")
    try:
        with requests.get(
            url,
            headers={"User-Agent": "vs-benchmark-rcsb-diagnostics/0.1", "Range": "bytes=0-0"},
            timeout=timeout,
            stream=True,
        ) as response:
            next(response.iter_content(chunk_size=1), b"")
            return MethodResult("Python requests.get", response.ok, str(response.status_code), "")
    except requests.RequestException as exc:
        return MethodResult("Python requests.get", False, "", str(exc))


def test_powershell(url: str, timeout: int) -> MethodResult:
    command = [
        "powershell",
        "-NoProfile",
        "-ExecutionPolicy",
        "Bypass",
        "-Command",
        (
            "$ErrorActionPreference='Stop'; "
            f"$ProgressPreference='SilentlyContinue'; "
            f"$r=Invoke-WebRequest -Uri {json.dumps(url)} -Method Head -UseBasicParsing -TimeoutSec {timeout}; "
            "[PSCustomObject]@{StatusCode=[int]$r.StatusCode; StatusDescription=$r.StatusDescription} | ConvertTo-Json -Compress"
        ),
    ]
    try:
        completed = subprocess.run(command, capture_output=True, text=True, timeout=timeout + 5)
        if completed.returncode == 0:
            status_code = parse_json_status(completed.stdout)
            return MethodResult("PowerShell Invoke-WebRequest", True, status_code, "")
        return MethodResult(
            "PowerShell Invoke-WebRequest",
            False,
            "",
            clean_subprocess_error(completed.stderr or completed.stdout),
        )
    except (subprocess.SubprocessError, OSError) as exc:
        return MethodResult("PowerShell Invoke-WebRequest", False, "", str(exc))


def test_curl(url: str, timeout: int, curl_path: str) -> MethodResult:
    command = [curl_path, "-I", "-L", "--max-time", str(timeout), "-sS", "-o", "NUL", "-w", "%{http_code}", url]
    try:
        completed = subprocess.run(command, capture_output=True, text=True, timeout=timeout + 5)
        status_code = completed.stdout.strip()
        if completed.returncode == 0 and status_code and status_code != "000":
            return MethodResult("curl -I", True, status_code, "")
        return MethodResult("curl -I", False, status_code if status_code != "000" else "", clean_subprocess_error(completed.stderr))
    except (subprocess.SubprocessError, OSError) as exc:
        return MethodResult("curl -I", False, "", str(exc))


def parse_json_status(text: str) -> str:
    try:
        parsed = json.loads(text)
        return str(parsed.get("StatusCode", ""))
    except json.JSONDecodeError:
        return ""


def clean_subprocess_error(text: str) -> str:
    return " ".join(text.replace("\r", " ").replace("\n", " ").split())


def find_curl() -> str | None:
    for name in ["curl.exe", "curl"]:
        path = shutil.which(name)
        if path:
            return path
    return None


def classify_issue(dns_success: bool, results: list[MethodResult]) -> str:
    if not dns_success:
        return "DNS resolution failed before any HTTP connection."
    successes = [result for result in results if result.success]
    if successes:
        failed_python = any(result.method.startswith("Python") and not result.success for result in results)
        non_python_success = any(
            result.success and (result.method.startswith("PowerShell") or result.method.startswith("curl")) for result in results
        )
        if failed_python and non_python_success:
            return "Python-specific networking/proxy/SSL issue is likely; PowerShell/curl fallback may help."
        return "At least one method can reach the URL; full downloader should be able to use a working method."
    errors = " | ".join(result.error for result in results if result.error).lower()
    if "actively refused" in errors or "10061" in errors or "connection refused" in errors:
        return "DNS resolves, but TCP/HTTPS connection is refused; proxy, firewall, VPN, or network policy is likely."
    if "ssl" in errors or "certificate" in errors:
        return "SSL/certificate verification problem is likely."
    if "proxy" in errors:
        return "Proxy configuration problem is likely."
    if "timed out" in errors or "timeout" in errors:
        return "Connection timeout; firewall, routing, or unstable network is likely."
    return "All methods failed; network/proxy/firewall restrictions are likely."


def write_report(path: Path, candidate_csv: Path, diagnostics: list[UrlDiagnostics], curl_path: str | None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tested_ids = sorted({item.pdb_id for item in diagnostics})
    all_methods_failed = all(not result.success for item in diagnostics for result in item.method_results)
    any_python_failed_non_python_success = any(
        any(result.method.startswith("Python") and not result.success for result in item.method_results)
        and any(
            result.success and (result.method.startswith("PowerShell") or result.method.startswith("curl"))
            for result in item.method_results
        )
        for item in diagnostics
    )

    lines = [
        "# RCSB download diagnostics report",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        f"Candidate CSV: `{candidate_csv}`",
        f"Tested PDB IDs: {', '.join(tested_ids)}",
        f"curl executable: `{curl_path}`" if curl_path else "curl executable: not found",
        "",
        "## Overall assessment",
        "",
    ]
    if all_methods_failed:
        lines.extend(
            [
                "- 当前本机无法直接连接 RCSB archive.",
                "- 建议更换网络、关闭/配置代理、使用浏览器手动下载、或在 AutoDL/Linux 环境下载.",
            ]
        )
    elif any_python_failed_non_python_success:
        lines.append("- Python download failed while PowerShell/curl succeeded for at least one URL; downloader fallback is recommended.")
    else:
        lines.append("- At least one diagnostic method succeeded; see per-URL details below.")

    for item in diagnostics:
        lines.extend(
            [
                "",
                f"## {item.pdb_id} {item.format_name}",
                "",
                f"- URL: `{item.url}`",
                f"- Host: `{item.host}`",
                f"- DNS: {'success' if item.dns_success else 'failed'}",
                f"- DNS addresses: `{item.dns_addresses}`" if item.dns_addresses else "- DNS addresses: none",
                f"- DNS error: `{item.dns_error}`" if item.dns_error else "- DNS error: none",
                f"- Likely issue: {item.likely_issue}",
                "",
                "| Method | Success | HTTP status | Error |",
                "| --- | --- | --- | --- |",
            ]
        )
        for result in item.method_results:
            lines.append(
                f"| {escape_md(result.method)} | {'yes' if result.success else 'no'} | "
                f"{escape_md(result.status_code or '')} | {escape_md(result.error or '')} |"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def escape_md(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", " ")


if __name__ == "__main__":
    raise SystemExit(main())
