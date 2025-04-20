import pandas as pd
import csv


RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
PREFIXES = {
    "domain": "d_",
    "phylum": "p_",
    "class": "c_",
    "order": "o_",
    "family": "f_",
    "genus": "g_",
    "species": "s_"
}


def fix_homonyms(input_path: str, output_path: str):
    columns = ["seqid", "taxname", "taxid"] + RANKS + ["sequence"]
    df = pd.read_csv(input_path, sep='\t', header=None, names=columns)

    for i in range(len(RANKS) - 1, 0, -1):
        current = RANKS[i]
        parent = RANKS[i - 1]
        homonym_counts = df.groupby(current)[parent].nunique()
        homonyms = homonym_counts[homonym_counts > 1].index
        print(f"Homonyms in {current}: {list(homonyms)}")

        mask = df[current].isin(homonyms)
        df.loc[mask, current] = df.loc[mask, parent] + '_' + df.loc[mask, current]

    df.to_csv(output_path, sep='\t', index=False, header=False)


def cleanup_name(name: str):
    name = name.strip()
    name = name.replace(" ", "_")
    if "sp." in name:
        name = "NA"
    return name


def build_lineage(row: dict):
    lineage = []
    last_good = None
    last_good_index = None

    for i, rank in enumerate(RANKS):
        name = row[rank].strip()
        name = cleanup_name(name)
        if name != "" and name != "NA":
            lineage.append(name)
            last_good = name
            last_good_index = i
        elif rank == "species":
            return None
        else:
            if last_good is None:
                lineage.append(PREFIXES[rank].rstrip("_"))
            else:
                missing_prefix = "_".join(PREFIXES[RANKS[j]].rstrip("_") for j in range(last_good_index + 1, i + 1))
                lineage.append(f"{missing_prefix}_{last_good}")
    return lineage


def tsv_to_filled_taxonomy_tsv(input_path: str, output_path: str):
    with open(input_path, newline="") as input_file, open(output_path, "w") as output_file:
        reader = csv.DictReader(input_file, delimiter="\t", fieldnames=["seqid", "taxname", "taxid", *RANKS, "sequence"])
        for row in reader:
            lineage = build_lineage(row)
            if lineage is not None:
                output_file.write("\t".join([row["seqid"], row["taxname"], row["taxid"]] + lineage + [row["sequence"]]) + "\n")


def generate_fasta_and_taxonomy(input_path: str, output_fasta_path: str, output_taxonomy_path: str):
    with open(input_path, newline="") as input_file, open(output_fasta_path, "w") as fasta_file, open(output_taxonomy_path, "w") as taxonomy_file:
        reader = csv.DictReader(input_file, delimiter="\t", fieldnames=[
            "seqid", "taxname", "taxid", *RANKS, "sequence"
        ])
        taxonomy_file.write("\t".join(["seqid"] + RANKS) + "\n")
        for row in reader:
            fasta_file.write(f">{row['seqid']}\troot;" + ';'.join([row[rank] for rank in RANKS]) + f"\n{row['sequence']}\n")
            taxonomy_file.write(f"{row['seqid']}\t" + "\t".join([row[rank] for rank in RANKS]) + "\n")
