from ednarefdb import PrimerSet, DatabaseBuilder, NucleotideDataset, ReferenceDatabase, NucleotideFile
import json
import logging
import os
import urllib.request


logging.basicConfig(level=logging.INFO)

WORKING_DIR = "./workdir"
TAXONOMY_FILES = ("names.dmp", "nodes.dmp", "nucl_gb.accession2taxid")


def primer_set_from_assay(code: str, **kwargs) -> PrimerSet:
    with urllib.request.urlopen("https://platform.ednaexpeditions.org/api/v1/assays/") as response:
        assays = json.load(response)
    assay = next(a for a in assays if a["code"] == code.upper())
    return PrimerSet(
        name=code.lower(),
        fwd=assay["forward_primer"],
        rev=assay["reverse_primer"],
        min_length=assay["min_amplicon_length"],
        max_length=assay["max_amplicon_length"],
        **kwargs,
    )


# Broadest available NCBI downloads per marker (12s_advanced includes full mitogenomes).
dataset_ncbi_coi = NucleotideDataset(
    name="ncbi_coi",
    files=[NucleotideFile(path="/Volumes/acasis/ncbi/coi_sequences.fasta", type="ncbi")],
)
dataset_ncbi_12s_advanced = NucleotideDataset(
    name="ncbi_12s_advanced",
    files=[NucleotideFile(path="/Volumes/acasis/ncbi/12s_advanced_sequences.fasta", type="ncbi")],
)
dataset_ncbi_16s_eukaryotic = NucleotideDataset(
    name="ncbi_16s_eukaryotic",
    files=[NucleotideFile(path="/Volumes/acasis/ncbi/16s_eukaryotic_sequences.fasta", type="ncbi")],
)

# Platform length bounds are unset for EE_* assays; use marker-typical ranges.
primer_set_ee_16s = primer_set_from_assay("EE_16S", mismatch=2, pga_percid=0.7)
primer_set_ee_16s.min_length = 50
primer_set_ee_16s.max_length = 300

primer_set_ee_coi = primer_set_from_assay("EE_COI", mismatch=4, pga_percid=0.7)
primer_set_ee_coi.min_length = 250
primer_set_ee_coi.max_length = 400

primer_set_ee_mifish = primer_set_from_assay("EE_MIFISH", mismatch=2, pga_percid=0.7)
primer_set_ee_mifish.min_length = 100
primer_set_ee_mifish.max_length = 250

primer_set_ee_mimammal = primer_set_from_assay("EE_MIMAMMAL", mismatch=2, pga_percid=0.7)
primer_set_ee_mimammal.min_length = 100
primer_set_ee_mimammal.max_length = 250

primer_set_ee_teleo = primer_set_from_assay("EE_TELEO", mismatch=2, pga_percid=0.7)
primer_set_ee_teleo.min_length = 40
primer_set_ee_teleo.max_length = 150

database_ee_16s = ReferenceDatabase(dataset_ncbi_16s_eukaryotic, primer_set_ee_16s)
database_ee_coi = ReferenceDatabase(dataset_ncbi_coi, primer_set_ee_coi)
database_ee_mifish = ReferenceDatabase(dataset_ncbi_12s_advanced, primer_set_ee_mifish)
database_ee_mimammal = ReferenceDatabase(dataset_ncbi_12s_advanced, primer_set_ee_mimammal)
database_ee_teleo = ReferenceDatabase(dataset_ncbi_12s_advanced, primer_set_ee_teleo)

os.makedirs(WORKING_DIR, exist_ok=True)
if not all(os.path.exists(os.path.join(WORKING_DIR, name)) for name in TAXONOMY_FILES):
    logging.info("NCBI taxonomy files missing from %s; downloading", WORKING_DIR)
    DatabaseBuilder(
        database_ee_16s, working_dir=WORKING_DIR, environment=None, dry_run=False
    ).ncbi_download_taxonomy()

for database in [
    database_ee_16s,
    # database_ee_coi,
    # database_ee_mifish,
    # database_ee_mimammal,
    # database_ee_teleo,
]:
    builder = DatabaseBuilder(database, working_dir=WORKING_DIR, environment=None, dry_run=False)
    builder.build()
    builder.compress_sintax()
