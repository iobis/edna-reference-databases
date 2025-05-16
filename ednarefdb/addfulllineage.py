import sys


def addfulllineage(input_path: str, fasta_path: str, output_path: str):

    f1 = open(input_path, "r").readlines()

    with open(output_path, "w") as output_file:

        hash = {}  # lineage map
        for line in f1[1:]:
            cols = line.strip().split("\t")
            lineage = ["Root"]
            for node in cols[1:]:
                if not node == "-":
                    lineage.append(node)
            ID = cols[0]
            lineage = ",".join(lineage)
            hash[ID] = lineage
        f2 = open(fasta_path, "r").readlines()
        for line in f2:
            if line[0] == ">":
                ID = line.strip().replace(">", "")
                try:
                    lineage = hash[ID]
                except KeyError:
                    output_file.write(ID, "not in taxonomy file" + "\n")
                    sys.exit()
                output_file.write(line.strip() + "\t" + lineage + "\n")
            else:
                output_file.write(line.strip() + "\n")
