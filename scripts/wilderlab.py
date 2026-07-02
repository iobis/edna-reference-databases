from ednarefdb import PrimerSet, DatabaseBuilder, NucleotideDataset, ReferenceDatabase, NucleotideFile
import json
import logging
import urllib.request


logging.basicConfig(level=logging.INFO)


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


dataset_ncbi_coi = NucleotideDataset(name="ncbi_coi", files=[
    NucleotideFile(path="/Volumes/acasis/ncbi/coi_sequences.fasta", type="ncbi")
])
dataset_ncbi_12s = NucleotideDataset(name="ncbi_12s", files=[
    NucleotideFile(path="/Volumes/acasis/ncbi/12s_sequences.fasta", type="ncbi")
])
dataset_ncbi_18s = NucleotideDataset(name="ncbi_18s", files=[
    NucleotideFile(path="/Volumes/acasis/ncbi/18s_sequences.fasta", type="ncbi")
])
# dataset_ncbi_euk_16s = NucleotideDataset(name="ncbi_euk_16s", files=[
#     NucleotideFile(path="/data/pieter/fasta/ncbi_euk_16s_50_50000.fasta", type="ncbi")
# ])

primer_set_ci = primer_set_from_assay("CI", mismatch=4, pga_percid=0.7)
primer_set_lx = primer_set_from_assay("LX", mismatch=4, pga_percid=0.7)
primer_set_rv = primer_set_from_assay("RV", mismatch=4, pga_percid=0.7)
primer_set_lm = primer_set_from_assay("LM", mismatch=0, pga_percid=0.7)
primer_set_bu = primer_set_from_assay("BU", mismatch=0, pga_percid=0.7)
primer_set_be = primer_set_from_assay("BE", mismatch=0, pga_percid=0.7)
primer_set_bx = primer_set_from_assay("BX", mismatch=0, pga_percid=0.7)
# primer_set_bu = PrimerSet(name="bu", fwd="TTGTACACACCGCCC", rev="CCTTCYGCAGGTTCACCTAC", mismatch=4, pga_percid=0.7, min_length=80, max_length=170)
# primer_set_wv = PrimerSet(name="wv", fwd="GACGAGAAGACCCTWTGGAGC", rev="CCRYGGTCGCCCCAAC", mismatch=4, pga_percid=0.7, min_length=35, max_length=200)
# primer_set_rv = PrimerSet(name="rv", fwd="TTAGATACCCCACTATGC", rev="TAGAACAGGCTCCTCTAG", mismatch=4, pga_percid=0.7, min_length=73, max_length=110)
# primer_set_hd = PrimerSet(name="hd", fwd="GGACGATAAGACCCTATAAA", rev="ACGCTGTTATCCCTAAAGT", mismatch=4, pga_percid=0.7, min_length=105, max_length=230)
# primer_set_bx = PrimerSet(name="bx", fwd="GCCAGTAGTCATATGCTTGTCT", rev="GCCTGCTGCCTTCCTT", mismatch=4, pga_percid=0.7, min_length=350, max_length=450)
# primer_set_gf = PrimerSet(name="gf", fwd="", rev="")
# primer_set_lm = PrimerSet(name="lm", fwd="CGTGCCAGCCACCGCG", rev="GGGTATCTAATCCYAGTTTG", mismatch=4, pga_percid=0.7, min_length=150, max_length=212)
# primer_set_um = PrimerSet(name="um", fwd="GGATTAGATACCCTGGTA", rev="CCGTCAATTCMTTTRAGTTT")

database_ci = ReferenceDatabase(dataset_ncbi_coi, primer_set_ci)
database_lx = ReferenceDatabase(dataset_ncbi_12s, primer_set_lx)
database_rv = ReferenceDatabase(dataset_ncbi_12s, primer_set_rv)
database_lm = ReferenceDatabase(dataset_ncbi_12s, primer_set_lm)
database_bu = ReferenceDatabase(dataset_ncbi_18s, primer_set_bu)
database_be = ReferenceDatabase(dataset_ncbi_18s, primer_set_be)
database_bx = ReferenceDatabase(dataset_ncbi_18s, primer_set_bx)
# database = ReferenceDatabase(dataset_ncbi_euk_16s, primer_set_wv)
# database = ReferenceDatabase(dataset_ncbi_12s, primer_set_rv)
# database = ReferenceDatabase(dataset_ncbi_euk_16s, primer_set_hd)
# database = ReferenceDatabase(dataset_ncbi_12s, primer_set_lm)
# database = ReferenceDatabase(dataset_ncbi_18s, primer_set_bx)

# builder.ncbi_download_taxonomy()

for database in [
    # database_ci,
    # database_rv,
    # database_lm,
    # database_bu,
    # database_be,
    # database_bx,
    database_lx,
]:
    builder = DatabaseBuilder(database, working_dir="./workdir", environment=None, dry_run=False)
    builder.build()
    builder.compress_sintax()
