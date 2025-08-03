from ednarefdb import PrimerSet, DatabaseBuilder, NucleotideDataset, ReferenceDatabase, NucleotideFile
import logging


logging.basicConfig(level=logging.INFO)


dataset_ncbi_coi = NucleotideDataset(name="ncbi_coi", files=[
    NucleotideFile(path="/data/pieter/fasta/ncbi_coi_50_50000.fasta", type="ncbi"),
    # NucleotideFile(path="/data/pieter/fasta/BOLD_Public.02-May-2025_coi.fasta", type="bold")
])
dataset_ncbi_euk_16s = NucleotideDataset(name="ncbi_euk_16s", files=[
    NucleotideFile(path="/data/pieter/fasta/ncbi_euk_16s_50_50000.fasta", type="ncbi")
])
dataset_ncbi_12s = NucleotideDataset(name="ncbi_12s", files=[
    NucleotideFile(path="/data/pieter/fasta/12s_50_50000.fasta", type="ncbi")
])
dataset_ncbi_18s = NucleotideDataset(name="ncbi_18s", files=[
    NucleotideFile(path="/data/pieter/fasta/18s_50_50000.fasta", type="ncbi")
])

primer_set_ci = PrimerSet(name="ci", fwd="DACWGGWTGAACWGTWTAYCCHCC", rev="GTTGTAATAAAATTAAYDGCYCCTARAATDGA", mismatch=6, pga_percid=0.7, min_length=76, max_length=76)
primer_set_wv = PrimerSet(name="wv", fwd="GACGAGAAGACCCTWTGGAGC", rev="CCRYGGTCGCCCCAAC", mismatch=4, pga_percid=0.7, min_length=35, max_length=200)
primer_set_rv = PrimerSet(name="rv", fwd="TTAGATACCCCACTATGC", rev="TAGAACAGGCTCCTCTAG", mismatch=4, pga_percid=0.7, min_length=73, max_length=110)
primer_set_hd = PrimerSet(name="hd", fwd="GGACGATAAGACCCTATAAA", rev="ACGCTGTTATCCCTAAAGT", mismatch=4, pga_percid=0.7, min_length=105, max_length=230)
primer_set_bu = PrimerSet(name="bu", fwd="TTGTACACACCGCCC", rev="CCTTCYGCAGGTTCACCTAC", mismatch=4, pga_percid=0.7, min_length=80, max_length=170)
primer_set_bx = PrimerSet(name="bx", fwd="GCCAGTAGTCATATGCTTGTCT", rev="GCCTGCTGCCTTCCTT", mismatch=4, pga_percid=0.7, min_length=350, max_length=450)
# primer_set_gf = PrimerSet(name="gf", fwd="", rev="")
primer_set_lm = PrimerSet(name="lm", fwd="CGTGCCAGCCACCGCG", rev="GGGTATCTAATCCYAGTTTG", mismatch=4, pga_percid=0.7, min_length=150, max_length=212)
# primer_set_um = PrimerSet(name="um", fwd="GGATTAGATACCCTGGTA", rev="CCGTCAATTCMTTTRAGTTT")

# database = ReferenceDatabase(dataset_ncbi_coi, primer_set_ci)
# database = ReferenceDatabase(dataset_ncbi_euk_16s, primer_set_wv)
# database = ReferenceDatabase(dataset_ncbi_12s, primer_set_rv)
# database = ReferenceDatabase(dataset_ncbi_euk_16s, primer_set_hd)
# database = ReferenceDatabase(dataset_ncbi_18s, primer_set_bu)
# database = ReferenceDatabase(dataset_ncbi_12s, primer_set_lm)
database = ReferenceDatabase(dataset_ncbi_18s, primer_set_bx)

builder = DatabaseBuilder(database, working_dir="/data/pieter/workdir", environment=None, dry_run=False)
# builder.ncbi_download_taxonomy()
builder.build()
