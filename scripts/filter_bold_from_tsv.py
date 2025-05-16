from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
import sys

input_file = "/data/pieter/fasta/BOLD_Public.02-May-2025.tsv"
output_file = "/data/pieter/fasta/BOLD_Public.02-May-2025_coi.fasta"

csv.field_size_limit(sys.maxsize)


def replace_non_ascii_with_underscore(s):
    return ''.join(c if ord(c) < 128 else '_' for c in s)


with open(output_file, "w") as f_out, open(input_file, "r") as f_in:

    reader = csv.DictReader(f_in, delimiter="\t")
    for i, line in enumerate(reader):
        if not i % 1000000:
            print(i)

        if line["marker_code"] == "COI-5P":
            if line["insdc_acs"] is not None and line["insdc_acs"] != "" and line["insdc_acs"] != "None":
                seqid = line["insdc_acs"].replace("-SUPPRESSED", "").replace("-WITHDRAWN", "").replace("ï»¿", "")
            else:
                seqid = line["processid"]

            header = "|".join([
                line["processid"],
                line["identification"],
                line["marker_code"],
                seqid
            ])

            record = SeqRecord(
                Seq(line["nuc"].replace("-", "")),
                id=replace_non_ascii_with_underscore(header),
                name=replace_non_ascii_with_underscore(seqid),
                description="",
            )

            SeqIO.write(record, f_out, "fasta")
