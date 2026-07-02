#!/usr/bin/env python3
"""
Download NCBI nucleotide sequences as FASTA using Biopython Entrez.

Uses Entrez history for large queries and parallel efetch batches.

Example:
    python scripts/download_ncbi_old.py --email you@example.com --preset coi
    python scripts/download_ncbi_old.py --email you@example.com --preset 18s
    python scripts/download_ncbi_old.py --email you@example.com --preset its
    python scripts/download_ncbi_old.py --email you@example.com --preset 12s
    python scripts/download_ncbi_old.py --email you@example.com --preset 16s_eukaryotic
    python scripts/download_ncbi_old.py --email you@example.com --preset 16s_prokaryotic
    python scripts/download_ncbi_old.py --email you@example.com --preset coi --workers 5 --api-key KEY
"""

from __future__ import annotations

import argparse
import json
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from Bio import Entrez

PRESET_QUERIES = {
    "coi": (
        "(COI[All Fields] OR CO1[All Fields] OR cytochrome oxidase subunit I[All Fields] "
        "OR cytochrome oxidase subunit 1[All Fields] OR cytochrome c oxidase subunit I[All Fields] "
        "OR cytochrome c oxidase subunit 1[All Fields] OR COX1[All Fields]) "
        'AND ("50"[SLEN] : "50000"[SLEN])'
    ),
    "18s": (
        '("Eukaryota"[Organism]) AND (18S[All Fields] OR "18S ribosomal RNA"[All Fields] '
        'OR "18S rRNA"[All Fields] OR "small subunit ribosomal RNA"[All Fields] OR SSU[Title]) '
        'AND ("50"[SLEN] : "50000"[SLEN]) NOT (mitochondrion[filter] OR chloroplast[filter])'
    ),
    "its": (
        '("Eukaryota"[Organism]) AND (ITS[All Fields] OR "internal transcribed spacer"[All Fields] '
        'OR ITS1[All Fields] OR ITS2[All Fields] OR "internal transcribed spacer 1"[All Fields] '
        'OR "internal transcribed spacer 2"[All Fields]) '
        'AND ("50"[SLEN] : "50000"[SLEN]) NOT (mitochondrion[filter] OR chloroplast[filter])'
    ),
    "12s": (
        'mitochondrion[filter] AND (12S[All Fields] OR "12S rRNA"[All Fields] '
        'OR "12S ribosomal RNA"[All Fields]) AND ribosomal[All Fields] '
        'AND ("50"[SLEN] : "50000"[SLEN])'
    ),
    "16s_eukaryotic": (
        '("Eukaryota"[Organism]) AND (16S[All Fields] OR "16S rRNA"[All Fields] '
        'OR "16S ribosomal RNA"[All Fields]) AND ribosomal[All Fields] '
        'AND mitochondrion[filter] AND ("50"[SLEN] : "50000"[SLEN])'
    ),
    "16s_prokaryotic": (
        '(Bacteria[Organism] OR Archaea[Organism]) AND (16S[All Fields] OR "16S rRNA"[All Fields] '
        'OR "16S ribosomal RNA"[All Fields]) AND ribosomal[All Fields] '
        'AND ("50"[SLEN] : "50000"[SLEN])'
    ),
}

DEFAULT_WORKERS = 3
DEFAULT_BATCH_SIZE = 500


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


def c(text: str, *codes: str) -> str:
    return f"{''.join(codes)}{text}{Colors.RESET}" if codes else text


def log(msg: str, level: str = "info") -> None:
    styles = {
        "info": (Colors.CYAN, "INFO"),
        "done": (Colors.GREEN, "DONE"),
        "warn": (Colors.YELLOW, "WARN"),
        "error": (Colors.RED, "ERROR"),
        "progress": (Colors.DIM, "...."),
        "header": (Colors.BOLD + Colors.BLUE, "===="),
    }
    code, tag = styles.get(level, (Colors.CYAN, "INFO"))
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"{c(f'[{timestamp}]', Colors.DIM)} {c(f'[{tag}]', code)} {msg}"
    with _print_lock:
        print(line, flush=True)


def human_duration(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.0f}s"
    if seconds < 3600:
        return f"{seconds / 60:.1f}m"
    return f"{seconds / 3600:.1f}h"


@dataclass
class BatchJob:
    batch_num: int
    retstart: int
    retmax: int


def progress_path(output: Path) -> Path:
    return output.with_name(output.name + ".progress")


def read_entrez_text(handle) -> str:
    data = handle.read()
    if isinstance(data, bytes):
        return data.decode("utf-8", errors="replace")
    return data


def is_ncbi_error_response(text: str) -> bool:
    stripped = text.lstrip()
    return not stripped or stripped.startswith("<?xml") or "<ERROR>" in stripped


class NcbiFastaDownloader:

    def __init__(
        self,
        email: str,
        api_key: str | None = None,
        workers: int = DEFAULT_WORKERS,
        batch_size: int = DEFAULT_BATCH_SIZE,
        max_retries: int = 10,
    ):
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        self.workers = workers
        self.batch_size = batch_size
        self.max_retries = max_retries
        self._write_lock = threading.Lock()
        self._progress_lock = threading.Lock()
        self._progress_file_lock = threading.Lock()
        self._sequences_written = 0
        self._batches_done = 0

    def esearch(self, query: str) -> tuple[int, str, str]:
        log(f"Searching NCBI nucleotide: {c(query, Colors.MAGENTA)}", level="header")
        handle = Entrez.esearch(db="nucleotide", term=query, usehistory="y", retmax=0)
        results = Entrez.read(handle)
        handle.close()
        count = int(results["Count"])
        log(f"Found {c(str(count), Colors.BOLD + Colors.BRIGHT_YELLOW)} records")
        return count, results["WebEnv"], results["QueryKey"]

    def load_completed_retstarts(self, output: Path) -> set[int]:
        path = progress_path(output)
        if not path.exists():
            return set()
        return set(json.loads(path.read_text()))

    def mark_batch_complete(self, output: Path, retstart: int) -> None:
        with self._progress_file_lock:
            completed = self.load_completed_retstarts(output)
            completed.add(retstart)
            progress_path(output).write_text(json.dumps(sorted(completed)))

    def efetch_fasta(self, webenv: str, query_key: str, retstart: int, retmax: int) -> str:
        kwargs = {
            "db": "nucleotide",
            "rettype": "fasta",
            "retmode": "text",
            "retstart": retstart,
            "retmax": retmax,
            "webenv": webenv,
            "query_key": query_key,
        }
        for attempt in range(self.max_retries):
            handle = None
            try:
                handle = Entrez.efetch(**kwargs)
                text = read_entrez_text(handle)
                if is_ncbi_error_response(text):
                    raise ValueError("NCBI returned an error document instead of FASTA")
                return text
            except Exception as e:
                if attempt == self.max_retries - 1:
                    raise
                wait = min(10 * (attempt + 1), 60)
                log(f"efetch retstart={retstart} failed ({e}), retry in {wait}s", level="warn")
                time.sleep(wait)
            finally:
                if handle is not None:
                    handle.close()
        return ""

    def count_sequences(self, fasta_text: str) -> int:
        return sum(1 for line in fasta_text.splitlines() if line.startswith(">"))

    def run_batch(self, job: BatchJob, webenv: str, query_key: str, output: Path) -> int:
        fasta_text = self.efetch_fasta(webenv, query_key, job.retstart, job.retmax)
        count = self.count_sequences(fasta_text)
        if fasta_text:
            with self._write_lock:
                with output.open("a") as out:
                    out.write(fasta_text if fasta_text.endswith("\n") else fasta_text + "\n")
        self.mark_batch_complete(output, job.retstart)
        return count

    def report_progress(
        self,
        batch_count: int,
        total_batches: int,
        batch_sequences: int,
        started: float,
    ) -> None:
        with self._progress_lock:
            self._sequences_written += batch_sequences
            self._batches_done += 1
            elapsed = time.monotonic() - started
            rate = self._sequences_written / elapsed if elapsed > 0 else 0
            remaining = total_batches - self._batches_done
            eta = (elapsed / self._batches_done) * remaining if self._batches_done else 0
            log(
                f"Batch {self._batches_done}/{total_batches}: "
                f"+{batch_sequences} seqs | "
                f"{c(str(self._sequences_written), Colors.BOLD)} total | "
                f"{rate:.0f} seq/s | "
                f"ETA {c(human_duration(eta), Colors.YELLOW)}",
                level="progress",
            )

    def download(
        self,
        query: str,
        output: Path,
        max_records: int | None = None,
        resume: bool = False,
    ) -> int:
        count, webenv, query_key = self.esearch(query)
        if max_records is not None:
            count = min(count, max_records)

        jobs = [
            BatchJob(batch_num=i + 1, retstart=start, retmax=min(self.batch_size, count - start))
            for i, start in enumerate(range(0, count, self.batch_size))
        ]
        total_batches = len(jobs)

        completed = self.load_completed_retstarts(output) if resume else set()
        if resume and completed:
            jobs = [job for job in jobs if job.retstart not in completed]
            self._sequences_written = sum(
                min(self.batch_size, count - retstart) for retstart in completed
            )
            self._batches_done = len(completed)
            log(
                f"Resuming: {len(completed)} batches already done, "
                f"{len(jobs)} remaining",
                level="header",
            )
        elif resume and output.exists() and output.stat().st_size > 0:
            log(
                "Resume requested but no progress file found; "
                "re-run without --resume to start over, or keep the partial file and "
                "fetch only the tail with a fresh query using --max-records",
                level="warn",
            )
        elif not resume:
            output.write_text("")
            progress = progress_path(output)
            if progress.exists():
                progress.unlink()

        if not jobs:
            log("Nothing left to download", level="done")
            return self._sequences_written

        log(
            f"Downloading {count} sequences with {c(str(self.workers), Colors.BOLD)} workers, "
            f"batch size {self.batch_size}",
            level="header",
        )
        started = time.monotonic()

        with ThreadPoolExecutor(max_workers=self.workers) as executor:
            futures = {
                executor.submit(self.run_batch, job, webenv, query_key, output): job
                for job in jobs
            }
            for future in as_completed(futures):
                job = futures[future]
                try:
                    batch_sequences = future.result()
                    self.report_progress(job.batch_num, total_batches, batch_sequences, started)
                except Exception as e:
                    log(f"Batch {job.batch_num} failed at retstart={job.retstart}: {e}", level="error")
                    raise

        elapsed = time.monotonic() - started
        log(
            f"Wrote {c(str(self._sequences_written), Colors.BOLD + Colors.GREEN)} sequences "
            f"to {output} in {human_duration(elapsed)}",
            level="done",
        )
        return self._sequences_written


def resolve_query(preset: str, query: str | None) -> str:
    if query is not None:
        return query
    return PRESET_QUERIES[preset]


def main():
    parser = argparse.ArgumentParser(description="Download NCBI nucleotide FASTA via Entrez")
    parser.add_argument("--email", required=True, help="Email for NCBI Entrez (required)")
    parser.add_argument("--api-key", help="NCBI API key (allows higher request rate)")
    parser.add_argument("--preset", choices=PRESET_QUERIES, default="coi", help="Preconfigured query")
    parser.add_argument("--query", help="Custom Entrez query (overrides --preset)")
    parser.add_argument("--output", type=Path, help="Output FASTA (default: {preset}_sequences.fasta)")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS, help="Parallel efetch workers")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    parser.add_argument("--max-records", type=int, help="Cap records downloaded")
    parser.add_argument("--resume", action="store_true", help="Resume a previous download")
    args = parser.parse_args()

    query = resolve_query(args.preset, args.query)
    output = args.output or Path(f"{args.preset}_sequences.fasta")

    downloader = NcbiFastaDownloader(
        email=args.email,
        api_key=args.api_key,
        workers=args.workers,
        batch_size=args.batch_size,
    )
    downloader.download(query, output, max_records=args.max_records, resume=args.resume)


if __name__ == "__main__":
    main()
