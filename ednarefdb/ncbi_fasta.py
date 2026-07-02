import logging
import os
import re
import subprocess
import tempfile
import glob
import shutil
from dataclasses import dataclass


@dataclass
class NcbiFastaBuildResult:
    fasta_path: str
    accessions_requested: int
    sequences_written: int


def batch(iterable, size):
    for i in range(0, len(iterable), size):
        yield iterable[i:i + size]


def read_accessions(path):
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def count_fasta_sequences(path):
    with open(path) as f:
        return sum(1 for line in f if line.startswith(">"))


def accession_from_header(header):
    token = header[1:].split()[0]
    if "|" in token:
        parts = token.split("|")
        for i, part in enumerate(parts[:-1]):
            if part in ("gb", "ref", "emb", "dbj", "pdb") and parts[i + 1]:
                return parts[i + 1]
    return token


def resolve_accession(header_acc, pending):
    if header_acc in pending:
        return header_acc
    if "." in header_acc:
        base = header_acc.rsplit(".", 1)[0]
        for acc in pending:
            if acc.rsplit(".", 1)[0] == base:
                return acc
    return None


def list_volume_tarballs(nt_dir):
    tarballs = []
    for path in glob.glob(os.path.join(nt_dir, "nt.*.tar.gz")):
        match = re.search(r"nt\.(\d+)\.tar\.gz$", os.path.basename(path))
        if match:
            tarballs.append((int(match.group(1)), path))
    return [path for _, path in sorted(tarballs)]


def validate_blast_db(nt_path: str, rolling_extract: bool = False):
    nt_dir = os.path.dirname(nt_path) or "."
    if glob.glob(f"{nt_path}.nhr") or glob.glob(f"{nt_path}.*.nhr"):
        return
    if rolling_extract and list_volume_tarballs(nt_dir):
        return
    raise FileNotFoundError(
        f"No BLAST database found at '{nt_path}' "
        f"(expected {nt_path}.nhr or {nt_path}.NNN.nhr). "
        f"If you have nt.*.tar.gz archives in {nt_dir}, use rolling_extract=True or "
        f"pass --rolling-extract to scripts/build_ncbi_fasta.py."
    )


class NcbiFastaBuilder:

    def __init__(
        self,
        nt_path: str,
        accessions_file: str,
        output_fasta: str,
        batch_size: int = 1000,
        rolling_extract: bool = False,
    ):
        self.nt_path = nt_path
        self.accessions_file = accessions_file
        self.output_fasta = output_fasta
        self.batch_size = batch_size
        self.rolling_extract = rolling_extract

    def build(self) -> NcbiFastaBuildResult:
        validate_blast_db(self.nt_path, self.rolling_extract)
        accessions = read_accessions(self.accessions_file)
        os.makedirs(os.path.dirname(self.output_fasta) or ".", exist_ok=True)

        with open(self.output_fasta, "w"):
            pass

        if self.rolling_extract:
            self.build_rolling(accessions)
        else:
            self.extract_accessions(accessions)

        sequences_written = count_fasta_sequences(self.output_fasta)

        logging.info(
            f"Wrote {sequences_written}/{len(accessions)} sequences to {self.output_fasta}"
        )

        return NcbiFastaBuildResult(
            fasta_path=self.output_fasta,
            accessions_requested=len(accessions),
            sequences_written=sequences_written,
        )

    def extract_accessions(self, accessions):
        pending = set(accessions)
        total_batches = (len(accessions) + self.batch_size - 1) // self.batch_size

        for i, accession_batch in enumerate(batch(accessions, self.batch_size), start=1):
            batch_found = self.run_blastdbcmd_batch(accession_batch, pending)
            pending -= batch_found
            logging.info(
                f"Batch {i}/{total_batches}: found {len(batch_found)}, "
                f"{len(pending)} still pending"
            )

    def build_rolling(self, accessions):
        nt_dir = os.path.dirname(self.nt_path) or "."
        tarballs = list_volume_tarballs(nt_dir)
        pending = set(accessions)

        for tarball in tarballs:
            if not pending:
                logging.info("All accessions found, stopping early")
                break

            volume_id = re.search(r"nt\.(\d+)\.tar\.gz$", os.path.basename(tarball)).group(1)
            logging.info(
                f"Volume {volume_id}: extracting ({len(pending)} accessions pending in memory)"
            )

            subprocess.run(["tar", "-xzf", tarball, "-C", nt_dir], check=True)
            try:
                volume_found = self.run_blastdbcmd_batches(list(pending), pending)
                pending -= volume_found
                logging.info(
                    f"Volume {volume_id}: found {len(volume_found)}, "
                    f"{len(pending)} still pending"
                )
            finally:
                self.remove_volume_files(nt_dir, volume_id)

    def run_blastdbcmd_batches(self, accessions, pending):
        found = set()
        total_batches = (len(accessions) + self.batch_size - 1) // self.batch_size

        for i, accession_batch in enumerate(batch(accessions, self.batch_size), start=1):
            batch_found = self.run_blastdbcmd_batch(accession_batch, pending)
            found |= batch_found
            if total_batches > 1:
                logging.info(f"  sub-batch {i}/{total_batches}: found {len(batch_found)}")

        return found

    def run_blastdbcmd_batch(self, accession_batch, pending):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as tmp:
            tmp.write("\n".join(accession_batch) + "\n")
            tmp_path = tmp.name

        batch_fasta = tempfile.NamedTemporaryFile(mode="w+", suffix=".fasta", delete=False)
        batch_fasta_path = batch_fasta.name
        batch_fasta.close()

        try:
            with open(batch_fasta_path, "w") as batch_out:
                subprocess.run(
                    [
                        "blastdbcmd", "-db", self.nt_path,
                        "-entry_batch", tmp_path, "-outfmt", "%f",
                    ],
                    check=True,
                    stdout=batch_out,
                )

            found = set()
            with open(batch_fasta_path) as batch_out:
                for line in batch_out:
                    if not line.startswith(">"):
                        continue
                    matched = resolve_accession(accession_from_header(line), pending)
                    if matched:
                        found.add(matched)

            if os.path.getsize(batch_fasta_path) > 0:
                with open(self.output_fasta, "a") as out, open(batch_fasta_path) as batch_in:
                    shutil.copyfileobj(batch_in, out)

            return found
        finally:
            os.remove(tmp_path)
            os.remove(batch_fasta_path)

    @staticmethod
    def remove_volume_files(nt_dir, volume_id):
        for path in glob.glob(os.path.join(nt_dir, f"nt.{volume_id}.n*")):
            os.remove(path)
