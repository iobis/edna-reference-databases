#!/usr/bin/env python3
"""Download all NCBI BLAST nt database files with fault recovery and progress reporting."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

BASE_URL = "https://ftp.ncbi.nlm.nih.gov/blast/db/"
NT_FILE_RE = re.compile(r"^nt(?:-[a-zA-Z0-9_.-]+|\.\d+\.tar\.gz(?:\.md5)?)$")
CHUNK_SIZE = 1024 * 1024  # 1 MiB
STATE_VERSION = 1
DEFAULT_WORKERS = 8
PROGRESS_INTERVAL = 2.0  # seconds between progress log lines per file
TOTAL_SUMMARY_INTERVAL = 60.0  # seconds between total-download summary lines


class Colors:
    RESET = "\033[0m"
    BOLD = "\033[1m"
    DIM = "\033[2m"
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    BRIGHT_YELLOW = "\033[93m"


_print_lock = threading.Lock()


def _c(text: str, *codes: str) -> str:
    if not codes:
        return text
    return f"{''.join(codes)}{text}{Colors.RESET}"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def human_bytes(n: int | float) -> str:
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    size = float(n)
    for unit in units:
        if size < 1024 or unit == units[-1]:
            return f"{size:.2f} {unit}" if unit != "B" else f"{int(size)} B"
        size /= 1024
    return f"{size:.2f} TiB"


def human_rate(bytes_per_sec: float) -> str:
    return f"{human_bytes(bytes_per_sec)}/s"


def log(
    msg: str,
    *,
    level: str = "info",
    also_print: bool = True,
) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    level_styles = {
        "info": (Colors.CYAN, "INFO"),
        "skip": (Colors.YELLOW, "SKIP"),
        "done": (Colors.GREEN, "DONE"),
        "error": (Colors.RED, "ERROR"),
        "warn": (Colors.YELLOW, "WARN"),
        "progress": (Colors.DIM, "...."),
        "header": (Colors.BOLD + Colors.BLUE, "===="),
        "total": (Colors.BOLD + Colors.BRIGHT_YELLOW, "TOTAL"),
    }
    code, tag = level_styles.get(level, (Colors.CYAN, "INFO"))
    if also_print:
        styled = f"{_c(f'[{timestamp}]', Colors.DIM)} {_c(f'[{tag}]', code)} {msg}"
        with _print_lock:
            print(styled, flush=True)


@dataclass
class FileEntry:
    name: str
    url: str
    expected_size: int | None = None


@dataclass
class DownloadState:
    version: int = STATE_VERSION
    started_at: str = field(default_factory=utc_now)
    updated_at: str = field(default_factory=utc_now)
    completed: dict[str, dict] = field(default_factory=dict)
    failed: dict[str, dict] = field(default_factory=dict)


class NtDownloader:
    def __init__(
        self,
        output_dir: Path,
        max_retries: int = 5,
        retry_delay: float = 5.0,
        verify_md5: bool = True,
        workers: int = DEFAULT_WORKERS,
    ) -> None:
        self.output_dir = output_dir
        self.max_retries = max_retries
        self.retry_delay = retry_delay
        self.verify_md5 = verify_md5
        self.workers = max(1, workers)
        self.state_path = output_dir / ".download_state.json"
        self.log_path = output_dir / "download.log"
        self.state = self._load_state()
        self._session_bytes = 0
        self._session_start = time.monotonic()
        self._lock = threading.Lock()
        self._progress_last: dict[str, float] = {}
        self._total_summary_last = 0.0

    def _load_state(self) -> DownloadState:
        if self.state_path.exists():
            data = json.loads(self.state_path.read_text())
            return DownloadState(
                version=data.get("version", STATE_VERSION),
                started_at=data.get("started_at", utc_now()),
                updated_at=data.get("updated_at", utc_now()),
                completed=data.get("completed", {}),
                failed=data.get("failed", {}),
            )
        return DownloadState()

    def _save_state(self) -> None:
        self.state.updated_at = utc_now()
        self.state_path.write_text(json.dumps(self.state.__dict__, indent=2) + "\n")

    def _append_log(self, msg: str) -> None:
        with self._lock:
            with self.log_path.open("a", encoding="utf-8") as fh:
                fh.write(f"[{utc_now()}] {msg}\n")

    def discover_files(self) -> list[FileEntry]:
        log("Fetching file list from NCBI FTP index...", level="header")
        req = urllib.request.Request(BASE_URL, headers={"User-Agent": "ncbi-nt-downloader/1.0"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            html = resp.read().decode("utf-8", errors="replace")

        entries: list[FileEntry] = []
        for match in re.finditer(
            r'<a href="([^"]+)">[^<]+</a>\s+[\d-]+\s+[\d:]+\s+([^\s<]+)',
            html,
        ):
            name, size_token = match.group(1), match.group(2)
            if not NT_FILE_RE.match(name):
                continue
            entries.append(
                FileEntry(
                    name=name,
                    url=BASE_URL + name,
                    expected_size=self._parse_size_token(size_token),
                )
            )

        entries.sort(key=lambda e: e.name)
        if not entries:
            raise RuntimeError("No nt files found in directory listing")
        log(f"Found {len(entries)} nt files to download", level="info")
        return entries

    @staticmethod
    def _parse_size_token(token: str) -> int | None:
        token = token.strip()
        if token == "-":
            return None
        multipliers = {"K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4}
        if token[-1] in multipliers:
            try:
                return int(float(token[:-1]) * multipliers[token[-1]])
            except ValueError:
                return None
        try:
            return int(token)
        except ValueError:
            return None

    def _remote_size(self, url: str) -> int | None:
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "ncbi-nt-downloader/1.0"})
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                length = resp.headers.get("Content-Length")
                return int(length) if length else None
        except urllib.error.HTTPError:
            return None

    def _local_md5(self, filename: str) -> str | None:
        md5_path = self.output_dir / f"{filename}.md5"
        if not md5_path.exists():
            return None
        line = md5_path.read_text().strip().split()[0]
        return line.lower()

    def _file_md5(self, path: Path) -> str:
        digest = hashlib.md5()
        with path.open("rb") as fh:
            while chunk := fh.read(CHUNK_SIZE):
                digest.update(chunk)
        return digest.hexdigest()

    def _cleanup_orphan_md5_files(self) -> int:
        removed = 0
        for md5_path in sorted(self.output_dir.glob("*.tar.gz.md5")):
            tar_path = self.output_dir / md5_path.name.removesuffix(".md5")
            if not tar_path.exists():
                md5_path.unlink()
                removed += 1
                self._append_log(f"removed orphan md5 {md5_path.name}")
        if removed:
            log(f"Removed {removed} orphan md5 file(s) with no matching archive", level="warn")
        return removed

    def _is_complete(self, entry: FileEntry) -> bool:
        dest = self.output_dir / entry.name
        if not dest.exists():
            return False

        size = dest.stat().st_size
        expected_size = entry.expected_size or self._remote_size(entry.url)
        if expected_size and size != expected_size:
            return False

        return True

    def _tar_is_current(self, tar: FileEntry) -> bool:
        dest = self.output_dir / tar.name
        if not dest.exists():
            return False

        if self.verify_md5:
            expected = self._local_md5(tar.name)
            if not expected:
                return False
            return self._file_md5(dest) == expected

        expected_size = tar.expected_size or self._remote_size(tar.url)
        if expected_size and dest.stat().st_size != expected_size:
            return False

        return True

    def _download_file(
        self,
        entry: FileEntry,
        file_idx: int,
        total_files: int,
        *,
        force: bool = False,
    ) -> None:
        dest = self.output_dir / entry.name
        part = dest.with_suffix(dest.suffix + ".part")
        expected_size = entry.expected_size or self._remote_size(entry.url)

        if not force and self._is_complete(entry):
            log(
                f"[{file_idx}/{total_files}] already complete: {entry.name}",
                level="skip",
            )
            with self._lock:
                self.state.completed[entry.name] = {
                    "size": dest.stat().st_size,
                    "completed_at": self.state.completed.get(entry.name, {}).get("completed_at", utc_now()),
                    "skipped": True,
                }
                self.state.failed.pop(entry.name, None)
                self._save_state()
            return

        if force and dest.exists():
            dest.unlink()
        if force and part.exists():
            part.unlink(missing_ok=True)

        attempt = 0
        while attempt < self.max_retries:
            attempt += 1
            resume_from = part.stat().st_size if part.exists() else 0
            if dest.exists() and resume_from == 0:
                dest.unlink()

            headers = {"User-Agent": "ncbi-nt-downloader/1.0"}
            if resume_from > 0:
                headers["Range"] = f"bytes={resume_from}-"

            req = urllib.request.Request(entry.url, headers=headers)
            started = time.monotonic()
            downloaded_this_attempt = 0

            try:
                with urllib.request.urlopen(req, timeout=300) as resp:
                    status = getattr(resp, "status", resp.getcode())
                    if resume_from > 0 and status not in (206, 200):
                        part.unlink(missing_ok=True)
                        raise urllib.error.HTTPError(
                            entry.url, status, "resume not supported", resp.headers, None
                        )
                    if status == 200 and resume_from > 0:
                        part.unlink(missing_ok=True)
                        resume_from = 0

                    mode = "ab" if resume_from > 0 else "wb"
                    with part.open(mode) as out:
                        while True:
                            chunk = resp.read(CHUNK_SIZE)
                            if not chunk:
                                break
                            out.write(chunk)
                            downloaded_this_attempt += len(chunk)
                            with self._lock:
                                self._session_bytes += len(chunk)
                            self._report_progress(
                                entry,
                                file_idx,
                                total_files,
                                part.stat().st_size,
                                expected_size,
                                started,
                            )

                part.rename(dest)
                if self.verify_md5 and entry.name.endswith(".tar.gz"):
                    expected = self._local_md5(entry.name)
                    if expected:
                        actual = self._file_md5(dest)
                        if actual != expected:
                            dest.unlink(missing_ok=True)
                            raise ValueError(f"MD5 mismatch for {entry.name}: {actual} != {expected}")

                elapsed = time.monotonic() - started
                rate = downloaded_this_attempt / elapsed if elapsed > 0 else 0
                log(
                    f"[{file_idx}/{total_files}] {entry.name} "
                    f"({human_bytes(dest.stat().st_size)}, {human_rate(rate)})",
                    level="done",
                )
                with self._lock:
                    self.state.completed[entry.name] = {
                        "size": dest.stat().st_size,
                        "completed_at": utc_now(),
                        "attempts": attempt,
                    }
                    self.state.failed.pop(entry.name, None)
                    self._save_state()
                self._append_log(f"completed {entry.name} after {attempt} attempt(s)")
                self._maybe_log_total_summary(force=entry.name.endswith(".tar.gz"))
                return

            except Exception as exc:  # noqa: BLE001 - retry loop needs broad catch
                delay = min(self.retry_delay * (2 ** (attempt - 1)), 300)
                msg = (
                    f"[{file_idx}/{total_files}] {entry.name} "
                    f"(attempt {attempt}/{self.max_retries}): {exc}"
                )
                log(msg, level="error")
                self._append_log(msg)
                with self._lock:
                    self.state.failed[entry.name] = {
                        "last_error": str(exc),
                        "attempts": attempt,
                        "last_attempt_at": utc_now(),
                    }
                    self._save_state()
                if attempt >= self.max_retries:
                    raise RuntimeError(f"Failed to download {entry.name} after {self.max_retries} attempts") from exc
                log(f"Retrying {entry.name} in {delay:.0f}s...", level="warn")
                time.sleep(delay)

    def _download_tar_bundle(
        self,
        tar: FileEntry,
        md5: FileEntry | None,
        md5_idx: int | None,
        tar_idx: int,
        total_files: int,
    ) -> None:
        if md5 is not None and md5_idx is not None:
            self._download_file(md5, md5_idx, total_files, force=True)
        tar_dest = self.output_dir / tar.name
        tar_part = tar_dest.with_suffix(tar_dest.suffix + ".part")
        tar_dest.unlink(missing_ok=True)
        tar_part.unlink(missing_ok=True)
        self._download_file(tar, tar_idx, total_files)

    def _total_archive_bytes(self) -> int:
        return sum(
            info.get("size", 0)
            for name, info in self.state.completed.items()
            if name.endswith(".tar.gz")
        )

    def _maybe_log_total_summary(self, *, force: bool = False) -> None:
        now = time.monotonic()
        with self._lock:
            if not force and now - self._total_summary_last < TOTAL_SUMMARY_INTERVAL:
                return
            self._total_summary_last = now
            session_bytes = self._session_bytes
            session_elapsed = max(now - self._session_start, 0.001)
            session_rate = session_bytes / session_elapsed
            archive_bytes = self._total_archive_bytes()
            archive_count = sum(
                1 for name in self.state.completed if name.endswith(".tar.gz")
            )

        log(
            f"Downloaded this session: {_c(human_bytes(session_bytes), Colors.BOLD + Colors.BRIGHT_YELLOW)} "
            f"@ {human_rate(session_rate)} | "
            f"Archives on disk: {_c(human_bytes(archive_bytes), Colors.BOLD + Colors.BRIGHT_YELLOW)} "
            f"({archive_count} complete)",
            level="total",
        )

    def _report_progress(
        self,
        entry: FileEntry,
        file_idx: int,
        total_files: int,
        current_size: int,
        expected_size: int | None,
        started: float,
    ) -> None:
        elapsed = max(time.monotonic() - started, 0.001)
        rate = current_size / elapsed if current_size else 0
        if expected_size and expected_size > 0:
            pct = 100.0 * current_size / expected_size
            remaining = expected_size - current_size
            eta = remaining / rate if rate > 0 else 0
            progress = (
                f"{human_bytes(current_size)} / {human_bytes(expected_size)} ({pct:.1f}%) "
                f"@ {human_rate(rate)} ETA {eta:.0f}s"
            )
        else:
            progress = f"{human_bytes(current_size)} @ {human_rate(rate)}"

        with self._lock:
            now = time.monotonic()
            last = self._progress_last.get(entry.name, 0.0)
            if now - last < PROGRESS_INTERVAL:
                return
            self._progress_last[entry.name] = now
            session_elapsed = max(time.monotonic() - self._session_start, 0.001)
            session_rate = self._session_bytes / session_elapsed
            completed_count = len(self.state.completed)

        log(
            f"[{file_idx}/{total_files}] {_c(entry.name, Colors.MAGENTA)}: {progress} | "
            f"session {human_bytes(self._session_bytes)} @ {human_rate(session_rate)} | "
            f"completed {completed_count}/{total_files}",
            level="progress",
        )
        self._maybe_log_total_summary()

    def run(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._cleanup_orphan_md5_files()
        files = self.discover_files()

        tar_files = [f for f in files if f.name.endswith(".tar.gz")]
        md5_by_tar = {
            f.name.removesuffix(".md5"): f for f in files if f.name.endswith(".md5")
        }
        meta_files = [f for f in files if not f.name.endswith((".tar.gz", ".md5"))]

        tars_to_download = [tar for tar in tar_files if not self._tar_is_current(tar)]
        tars_current = [tar for tar in tar_files if self._tar_is_current(tar)]
        total = len(meta_files) + len(tars_to_download) * 2

        for tar in tars_current:
            dest = self.output_dir / tar.name
            self.state.completed[tar.name] = {
                "size": dest.stat().st_size,
                "completed_at": self.state.completed.get(tar.name, {}).get("completed_at", utc_now()),
                "skipped": True,
            }
            self.state.failed.pop(tar.name, None)
            md5_name = f"{tar.name}.md5"
            if (self.output_dir / md5_name).exists():
                self.state.completed[md5_name] = {
                    "size": (self.output_dir / md5_name).stat().st_size,
                    "completed_at": self.state.completed.get(md5_name, {}).get("completed_at", utc_now()),
                    "skipped": True,
                }
                self.state.failed.pop(md5_name, None)
        self._save_state()

        total_bytes = sum(f.expected_size or 0 for f in tars_to_download)
        log(f"Output directory: {self.output_dir}", level="header")
        log(
            f"Archives: {len(tar_files)} total, {len(tars_to_download)} need download, "
            f"{len(tars_current)} already verified"
        )
        log(
            f"Planned transfers: {total} ({len(meta_files)} metadata, "
            f"{len(tars_to_download)} md5, {len(tars_to_download)} archives)"
        )
        if total_bytes:
            log(f"Estimated download size: {human_bytes(total_bytes)}")
        log(
            f"Starting downloads with {self.workers} worker(s) (Ctrl+C safe — rerun to resume)",
            level="header",
        )
        self._maybe_log_total_summary(force=True)

        jobs: list[tuple] = []
        idx = 0
        for entry in meta_files:
            idx += 1
            jobs.append(("file", entry, idx, total, False))

        for tar in tars_to_download:
            md5 = md5_by_tar.get(tar.name)
            md5_idx = None
            if md5 is not None:
                idx += 1
                md5_idx = idx
            idx += 1
            tar_idx = idx
            jobs.append(("tar_bundle", tar, md5, md5_idx, tar_idx, total))

        try:
            with ThreadPoolExecutor(max_workers=self.workers) as executor:
                futures = []
                for job in jobs:
                    if job[0] == "file":
                        _, entry, file_idx, job_total, force = job
                        futures.append(
                            executor.submit(self._download_file, entry, file_idx, job_total, force=force)
                        )
                    else:
                        _, tar, md5, md5_idx, tar_idx, job_total = job
                        futures.append(
                            executor.submit(
                                self._download_tar_bundle,
                                tar,
                                md5,
                                md5_idx,
                                tar_idx,
                                job_total,
                            )
                        )
                for future in as_completed(futures):
                    future.result()
        finally:
            with self._lock:
                self._save_state()

        failed = [name for name in self.state.failed if name not in self.state.completed]
        completed = len(self.state.completed)
        log(f"Finished: {completed}/{total} files complete", level="header")
        if failed:
            log(f"Failed files ({len(failed)}): {', '.join(failed)}", level="error")
            raise SystemExit(1)


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "nt",
        help="Directory to store downloaded files (default: ./nt)",
    )
    parser.add_argument(
        "--max-retries",
        type=int,
        default=5,
        help="Maximum retry attempts per file (default: 5)",
    )
    parser.add_argument(
        "--retry-delay",
        type=float,
        default=5.0,
        help="Initial retry delay in seconds, doubles each attempt (default: 5)",
    )
    parser.add_argument(
        "--no-md5-verify",
        action="store_true",
        help="Skip MD5 verification for .tar.gz files",
    )
    parser.add_argument(
        "-j",
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help=f"Number of parallel download threads (default: {DEFAULT_WORKERS})",
    )
    return parser.parse_args(list(argv) if argv is not None else None)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    downloader = NtDownloader(
        output_dir=args.output_dir.resolve(),
        max_retries=args.max_retries,
        retry_delay=args.retry_delay,
        verify_md5=not args.no_md5_verify,
        workers=args.workers,
    )
    downloader.run()


if __name__ == "__main__":
    main()
