from ednarefdb import PrimerSet, DatabaseBuilder, NucleotideDataset, ReferenceDatabase, NucleotideFile
import logging


logging.basicConfig(level=logging.INFO)


dataset_bold_ncbi_coi = NucleotideDataset(name="bold_ncbi_coi", files=[
    NucleotideFile(path="/data/pieter/fasta/ncbi_coi_50_50000.fasta", type="ncbi"),
    NucleotideFile(path="/data/pieter/fasta/BOLD_Public.02-May-2025_coi.fasta", type="bold")
])
primer_set_coi = PrimerSet(name="coi", fwd="GGWACWGGWTGAACWGTWTAYCCYCC", rev="TAIACYTCIGGRTGICCRAARAAYCA", min_length=200, max_length=50000, max_n=1)
database = ReferenceDatabase(dataset_bold_ncbi_coi, primer_set_coi)
builder = DatabaseBuilder(database, working_dir="/data/pieter/workdir", environment=None)

# builder.import_nucleotide()
# builder.pcr()
builder.pga(percid=0.8, coverage=0.9)
# builder.dereplicate()
# builder.filter()
# builder.export()
# builder.prepare_train()
# builder.train()
# builder.cleanup()






######### test NCBI download

# dataset_abra_coi = NucleotideDataset(name="abra_coi", type="ncbi", query='Abra[Organism] AND COI[All Fields] AND ("50"[SLEN] : "50000"[SLEN])')
# primer_set_coi = PrimerSet(name="coi", fwd="GGWACWGGWTGAACWGTWTAYCCYCC", rev="TAIACYTCIGGRTGICCRAARAAYCA", min_length=200, max_length=50000, max_n=1)
# database = ReferenceDatabase(dataset_abra_coi, primer_set_coi)
# builder = DatabaseBuilder(database, working_dir="/data/pieter/workdir", environment=None)
# builder.ncbi_download_taxonomy()
# builder.ncbi_download_nucleotide(dataset=dataset_abra_coi)
# builder.import_nucleotide()
