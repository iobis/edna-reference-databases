from ednarefdb import PrimerSet, DatabaseBuilder, NucleotideDataset, ReferenceDatabase
import logging


logging.basicConfig(level=logging.INFO)


primer_set_coi = PrimerSet(name="coi", fwd="GGWACWGGWTGAACWGTWTAYCCYCC", rev="TAIACYTCIGGRTGICCRAARAAYCA")
dataset_coi = NucleotideDataset(name="ncbi_coi_50_50000", path="/Volumes/acasis/fasta/ncbi_coi_50_50000.fasta", query='COI[All Fields] OR CO1[All Fields] OR cytochrome oxidase subunit I[All Fields] OR cytochrome oxidase subunit 1[All Fields] OR cytochrome c oxidase subunit I[All Fields] OR cytochrome c oxidase subunit 1[All Fields] OR COX1[All Fields] AND ("50"[SLEN] : "50000"[SLEN])')
dataset_ribo = NucleotideDataset(name="ncbi_ribosomal_50_50000", path="/Volumes/acasis/fasta/ncbi_ribosomal_50_50000.fasta", query='(12S[All Fields] OR 16S[All Fields]) AND ribosomal[All Fields] AND ("50"[SLEN] : "50000"[SLEN])')

datasets = [
    dataset_coi,
    # dataset_ribo
]

databases = [
    ReferenceDatabase(dataset_coi, primer_set_coi),
]

ref_db = DatabaseBuilder(working_dir="/Volumes/acasis/refdb_workdir", environment=None)
# ref_db.ncbi_download_taxonomy()

for database in databases:
    # ref_db.import_nucleotide(dataset=database.dataset)
    # ref_db.pcr(dataset=database.dataset, primer_set=database.primer_set)
    ref_db.pga(dataset=database.dataset, primer_set=database.primer_set, percid=0.8, coverage=0.8)
    # ref_db.sequence_cleanup(dataset=database.dataset, primer_set=database.primer_set)
    # ref_db.train(dataset=database.dataset, primer_set=database.primer_set)
    # ref_db.cleanup(dataset=database.dataset, primer_set=database.primer_set)
