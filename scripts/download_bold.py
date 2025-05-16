from ednarefdb import PrimerSet, DatabaseBuilder, NucleotideDataset, ReferenceDatabase, NucleotideFile
import logging


logging.basicConfig(level=logging.INFO)


ref_db = DatabaseBuilder(working_dir="/data/pieter/workdir", environment=None)

phyla = [
    "Acanthocephala",
    "Annelida",
    "Arthropoda",
    "Brachiopoda",
    "Bryozoa",
    "Chaetognatha",
    "Chordata",
    "Cnidaria",
    "Ctenophora",
    "Cycliophora",
    "Echinodermata",
    "Entoprocta",
    "Gastrotricha",
    "Gnathostomulida",
    "Hemichordata",
    "Kinorhyncha",
    "Mollusca",
    "Nematoda",
    "Nematomorpha",
    "Nemertea",
    "Onychophora",
    "Phoronida",
    "Placozoa",
    "Platyhelminthes",
    "Porifera",
    "Priapulida",
    "Rhombozoa",
    "Rotifera",
    "Tardigrada",
    "Xenacoelomorpha",
    "Bryophyta",
    "Chlorophyta",
    "Ascomycota",
    "Basidiomycota",
    "Chytridiomycota",
    "Glomeromycota",
    "Myxomycota",
    "Zygomycota",
    "Chlorarachniophyta",
    "Ciliophora",
    "Heterokontophyta",
    "Pyrrophycophyta",
    "Rhodophyta"
]

for phylum in phyla:
    print(f"Downloading phylum {phylum}")
    ref_db.run_command(f"""
        crabs --download-bold \
        --taxon '{phylum}' \
        --output temp_bold_{phylum}.fasta
    """)
