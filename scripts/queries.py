from ednarefdb import PrimerSet, DatabaseBuilder, NucleotideDataset, ReferenceDatabase
import logging


logging.basicConfig(level=logging.INFO)


primer_set_coi = PrimerSet(name="coi", fwd="GGWACWGGWTGAACWGTWTAYCCYCC", rev="TAIACYTCIGGRTGICCRAARAAYCA")
# primer_set_mifish = PrimerSet(name="mifish", fwd="GTYGGTAAAWCTCGTGCCAGC", rev="CATAGTGGGGTATCTAATCCYAGTTTG")

dataset_coi = NucleotideDataset(name="ncbi_coi_50_50000", path="/data/pieter/fasta/ncbi_coi_50_50000.fasta", query='COI[All Fields] OR CO1[All Fields] OR cytochrome oxidase subunit I[All Fields] OR cytochrome oxidase subunit 1[All Fields] OR cytochrome c oxidase subunit I[All Fields] OR cytochrome c oxidase subunit 1[All Fields] OR COX1[All Fields] AND ("50"[SLEN] : "50000"[SLEN])')
# dataset_12s = NucleotideDataset(name="ncbi_12s_50_50000", path="/data/pieter/fasta/ncbi_12s_50_50000.fasta", query='(12S[All Fields] OR 16S[All Fields]) AND ribosomal[All Fields] AND ("50"[SLEN] : "50000"[SLEN])')

datasets = [
    dataset_coi,
    # dataset_12s
]

databases = [
    ReferenceDatabase(dataset_coi, primer_set_coi),
    # ReferenceDatabase(dataset_12s, primer_set_mifish)
]

ref_db = DatabaseBuilder(working_dir="/data/pieter/refdb_workdir", environment=None)
# ref_db.ncbi_download_taxonomy()

for database in databases:
    # ref_db.import_nucleotide(dataset=database.dataset)
    # ref_db.pcr(dataset=database.dataset, primer_set=database.primer_set)
    # ref_db.pga(dataset=database.dataset, primer_set=database.primer_set, percid=0.8, coverage=0.8)
    # ref_db.dereplicate(dataset=database.dataset, primer_set=database.primer_set)
    ref_db.export(dataset=database.dataset, primer_set=database.primer_set)
    # ref_db.prepare_train(dataset=database.dataset, primer_set=database.primer_set)
    # ref_db.train(dataset=database.dataset, primer_set=database.primer_set)
    # ref_db.cleanup(dataset=database.dataset, primer_set=database.primer_set)
