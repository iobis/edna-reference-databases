from ednarefdb import PrimerSet, DatabaseBuilder, NucleotideDataset, ReferenceDatabase
import logging


logging.basicConfig(level=logging.INFO)


primer_set_coi = PrimerSet(name="coi", fwd="GGWACWGGWTGAACWGTWTAYCCYCC", rev="TAIACYTCIGGRTGICCRAARAAYCA")
dataset_coi = NucleotideDataset(name="coi_ncbi_1_50000", query='COI[All Fields] OR CO1[All Fields] OR cytochrome oxidase subunit I[All Fields] OR cytochrome oxidase subunit 1[All Fields] OR cytochrome c oxidase subunit I[All Fields] OR cytochrome c oxidase subunit 1[All Fields] OR COX1[All Fields] AND ("50"[SLEN] : "50000"[SLEN])')
dataset_ribo = NucleotideDataset(name="ncbi_ribo_1_50000", query='(12S[All Fields] OR 16S[All Fields]) AND ribosomal[All Fields] AND ("1"[SLEN] : "50000"[SLEN])')

datasets = [
    dataset_coi,
    # dataset_ribo
]

databases = [
    ReferenceDatabase(dataset_coi, primer_set_coi),
]


ref_db = DatabaseBuilder(working_dir="./output", environment=None)
# ref_db.ncbi_download_taxonomy()

for dataset in datasets:
    ref_db.ncbi_download_nucleotide(dataset)

# for database in databases:
    # ref_db.pcr(dataset=database.dataset, primer_set=database.primer_set)
    # ref_db.pga(dataset=database.dataset, primer_set=database.primer_set, percid=0.8, coverage=0.8)
    # ref_db.assign_taxonomy(dataset=database.dataset, primer_set=database.primer_set)
    # ref_db.sequence_cleanup(dataset=database.dataset, primer_set=database.primer_set)
    # ref_db.train(dataset=database.dataset, primer_set=database.primer_set)
    # ref_db.cleanup(dataset=database.dataset, primer_set=database.primer_set)
