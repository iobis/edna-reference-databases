from Bio import SeqIO

input_file = "/data/pieter/fasta/BOLD_Public.02-May-2025.fasta"
output_file = "/data/pieter/fasta/BOLD_Public.02-May-2025_coi.fasta"

with open(output_file, "w") as out_f:
    for i, record in enumerate(SeqIO.parse(input_file, "fasta")):
        if not i % 10000:
            print(i)
        parts = record.description.split("|")
        seqid = parts[0]
        marker = parts[1]
        location = parts[2]
        taxonomy = parts[3]

        lineage = [taxon for taxon in taxonomy.split(",") if taxon != "None"]
        if len(lineage) == 0:
            continue
        taxon = lineage[-1]

        # if marker == "COI-5P":
        if marker == "COI-5P":
            header = "|".join([seqid, taxon, marker, seqid])
            record.id = header
            record.description = ""
            SeqIO.write(record, out_f, "fasta")
