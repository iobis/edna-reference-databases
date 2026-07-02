import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ednarefdb import NcbiFastaBuilder


def main():
    parser = argparse.ArgumentParser(description="Build a FASTA from local nt and an accessions file")
    parser.add_argument("--nt", required=True, help="BLAST db prefix, e.g. /data/pieter/data/nt/nt")
    parser.add_argument("--accessions", required=True, help="Text file with one accession per line")
    parser.add_argument("--output", required=True, help="Output FASTA path")
    parser.add_argument("--batch-size", type=int, default=1000)
    parser.add_argument(
        "--rolling-extract",
        action="store_true",
        help="Extract one nt.*.tar.gz volume at a time; keep pending accessions in memory",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    result = NcbiFastaBuilder(
        nt_path=args.nt,
        accessions_file=args.accessions,
        output_fasta=args.output,
        batch_size=args.batch_size,
        rolling_extract=args.rolling_extract,
    ).build()

    print(
        f"Done: {result.sequences_written}/{result.accessions_requested} sequences "
        f"written to {result.fasta_path}"
    )


if __name__ == "__main__":
    main()
