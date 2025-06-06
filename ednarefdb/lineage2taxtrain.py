def lineage2taxtrain(input_path: str, output_path: str):
    """
    Used to convert a taxonomy in tab-delimited file containing the taxonomic hierarchical structure to RDP Classifier taxonomy training file
    Approach:each taxon is uniquely identified by the combination of its tax id and depth from the root rank, its attributes comprise: name, parent taxid, and level of depth from the root rank. 
    """

    f = open(input_path, "r").readlines()

    with open(output_path, "w") as output_file:

        header = f[0]
        cols = header.strip().split('\t')[1:]
        hash = {}  # taxon name-id map
        ranks = {}  # column number-rank map
        # lineages = []  # list of unique lineages

        hash = {"Root": 0}  # initiate root rank taxon id map
        for i in range(len(cols)):
            name = cols[i]
            ranks[i] = name
        root = ["0", "Root", "-1", "0", "rootrank"]  # root rank info
        output_file.write("*".join(root) + "\n")
        ID = 0  # taxon id
        for line in f[1:]:
            cols = line.strip().split('\t')[1:]
            # if not cols in lineages:#unique lineage
            # 	lineages.append(cols)
            for i in range(len(cols)):  # iterate each column
                # name = string.join(cols[:i + 1], ',')
                name = []
                for node in cols[:i + 1]:
                    if not node == "-":
                        name.append(node)
                pName = ",".join(name[:-1])
                # if not name in lineages:
                # 	lineages.append(name)
                depth = len(name)
                name = ",".join(name)
                if name in hash.keys():
                    continue
                rank = ranks[i]
                # level = len(name.split(','))
                # pName = string.join(cols[:i], ',')#parent name
                if i == 0:
                    pName = "Root"
                pID = hash[pName]  # parent taxid
                ID += 1
                hash[name] = ID  # add name-id to the map
                # out = ['%s'%ID, name, '%s'%pID, '%s'%depth, rank]
                out = ["%s" % ID, name.split(",")[-1], '%s' % pID, '%s' % depth, rank]
                output_file.write("*".join(out) + "\n")
