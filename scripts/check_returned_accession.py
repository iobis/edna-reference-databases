import subprocess

def get_blast_accession(accession, db="/data/pieter/data/nt/nt"):
    try:
        result = subprocess.run(
            ["blastdbcmd", "-db", db, "-entry", accession],
            capture_output=True, text=True, check=True
        )
        first_line = result.stdout.splitlines()[0]
        if first_line.startswith(">"):
            parsed_acc = first_line.split()[0][1:]
            return parsed_acc
        else:
            return None
    except subprocess.CalledProcessError as e:
        print(f"Error fetching {accession}: {e}")
        return None

with open("/data/pieter/data/coi_accessions.txt") as f:
    for line in f:
        acc = line.strip()
        if not acc:
            continue
        returned_acc = get_blast_accession(acc)
        if returned_acc is None:
            print(f"{acc}: Not found")
        elif returned_acc != acc:
            print(f"{acc}: Mismatch (got {returned_acc})")
        else:
            print(f"{acc}: Match")
