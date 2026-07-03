"""
Download nucleotide accessions from NCBI using Biopython Entrez.

Writes one accession per line. Optionally fetches matching sequences as FASTA.

Example:
    python scripts/download_ncbi_accessions.py --email you@example.com --preset coi
    python scripts/download_ncbi_accessions.py --email you@example.com --preset 12s
    python scripts/download_ncbi_accessions.py --email you@example.com --preset 16s_eukaryotic
    python scripts/download_ncbi_accessions.py --email you@example.com --preset 16s_prokaryotic
    python scripts/download_ncbi_accessions.py --email you@example.com --preset 28s
    python scripts/download_ncbi_accessions.py --email you@example.com --preset coi --fetch-sequences
    python scripts/download_ncbi_accessions.py --email you@example.com --query '16S[All Fields] AND ...'
"""

import argparse
import time
from pathlib import Path

from Bio import Entrez, SeqIO

PRESET_QUERIES = {
    "coi": (
        '(COI[All Fields] OR CO1[All Fields] OR cytochrome oxidase subunit I[All Fields] '
        'OR cytochrome oxidase subunit 1[All Fields] OR cytochrome c oxidase subunit I[All Fields] '
        'OR cytochrome c oxidase subunit 1[All Fields] OR COX1[All Fields]) '
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
    "28s": (
        '("Eukaryota"[Organism]) AND (28S[All Fields] OR "28S rRNA"[All Fields] '
        'OR "28S ribosomal RNA"[All Fields] OR "large subunit ribosomal RNA"[All Fields] OR LSU[Title]) '
        'AND ribosomal[All Fields] AND ("50"[SLEN] : "50000"[SLEN]) '
        "NOT (mitochondrion[filter] OR chloroplast[filter])"
    ),
}


def esearch_with_history(query):
    handle = Entrez.esearch(db="nucleotide", term=query, usehistory="y", retmax=0)
    results = Entrez.read(handle)
    handle.close()
    return int(results["Count"]), results["WebEnv"], results["QueryKey"]


def efetch_with_retry(efetch_kwargs, parse_fn, max_retries=10):
    for attempt in range(max_retries):
        handle = None
        try:
            handle = Entrez.efetch(**efetch_kwargs)
            return parse_fn(handle)
        except Exception as e:
            if attempt == max_retries - 1:
                raise
            wait = min(10 * (attempt + 1), 60)
            print(f"efetch failed ({e}), retrying in {wait}s...")
            time.sleep(wait)
        finally:
            if handle is not None:
                handle.close()


def fetch_accessions_batch(webenv, query_key, retstart, retmax):
    def parse(handle):
        return [line.strip() for line in handle if line.strip()]

    return efetch_with_retry(
        {
            "db": "nucleotide",
            "rettype": "acc",
            "retmode": "text",
            "retstart": retstart,
            "retmax": retmax,
            "webenv": webenv,
            "query_key": query_key,
        },
        parse,
    )


def fetch_fasta_batch(webenv, query_key, retstart, retmax):
    def parse(handle):
        return list(SeqIO.parse(handle, "fasta"))

    return efetch_with_retry(
        {
            "db": "nucleotide",
            "rettype": "fasta",
            "retmode": "text",
            "retstart": retstart,
            "retmax": retmax,
            "webenv": webenv,
            "query_key": query_key,
        },
        parse,
    )


def resolve_query(preset, query):
    if query is not None:
        return query
    return PRESET_QUERIES[preset]


def main():
    parser = argparse.ArgumentParser(description="Download NCBI nucleotide accessions via Entrez")
    parser.add_argument("--email", required=True, help="Email for NCBI Entrez (required by NCBI)")
    parser.add_argument("--preset", choices=PRESET_QUERIES, default="coi", help="Preconfigured Entrez query")
    parser.add_argument("--query", help="Custom Entrez query (overrides --preset)")
    parser.add_argument("--accessions-out", type=Path, help="Output accessions file (default: {preset}_accessions.txt)")
    parser.add_argument("--fasta-out", type=Path, help="Output FASTA file (default: {preset}_sequences.fasta)")
    parser.add_argument("--fetch-sequences", action="store_true", help="Also download sequences as FASTA")
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument("--max-records", type=int, default=None, help="Cap records fetched (default: all)")
    args = parser.parse_args()

    query = resolve_query(args.preset, args.query)
    accessions_out = args.accessions_out or Path(f"{args.preset}_accessions.txt")
    fasta_out = args.fasta_out or Path(f"{args.preset}_sequences.fasta")

    Entrez.email = args.email

    print(f"Searching NCBI nucleotide: {query}")
    count, webenv, query_key = esearch_with_history(query)
    if args.max_records is not None:
        count = min(count, args.max_records)
    total_batches = (count + args.batch_size - 1) // args.batch_size
    print(f"Found {count} records.")

    accessions_written = 0
    with accessions_out.open("w") as out:
        for batch_num, retstart in enumerate(range(0, count, args.batch_size), start=1):
            retmax = min(args.batch_size, count - retstart)
            accessions = fetch_accessions_batch(webenv, query_key, retstart, retmax)
            out.write("\n".join(accessions) + "\n")
            accessions_written += len(accessions)
            print(f"Accessions batch {batch_num}/{total_batches}: {accessions_written} total")

    print(f"Wrote {accessions_written} accessions to {accessions_out}")

    if not args.fetch_sequences:
        return

    sequences_written = 0
    with fasta_out.open("w") as out:
        for batch_num, retstart in enumerate(range(0, count, args.batch_size), start=1):
            retmax = min(args.batch_size, count - retstart)
            records = fetch_fasta_batch(webenv, query_key, retstart, retmax)
            sequences_written += SeqIO.write(records, out, "fasta")
            print(f"FASTA batch {batch_num}/{total_batches}: {sequences_written} sequences total")

    print(f"Wrote {sequences_written} sequences to {fasta_out}")


if __name__ == "__main__":
    main()
