from ednarefdb import PrimerSet, DatabaseBuilder, NucleotideDataset, ReferenceDatabase, NucleotideFile
import logging


logging.basicConfig(level=logging.INFO)


dataset_ncbi_coi = NucleotideDataset(name="ncbi_coi", files=[NucleotideFile(path="/data/pieter/fasta/ncbi_coi_50_50000.fasta", type="ncbi")])

primer_set_ci = PrimerSet(name="ci", fwd="DACWGGWTGAACWGTWTAYCCHCC", rev="GTTGTAATAAAATTAAYDGCYCCTARAATDGA", mismatch=8, pga_percid=0.7, min_length=76, max_length=76)
primer_set_bu = PrimerSet(name="bu", fwd="TTGTACACACCGCCC", rev="CCTTCYGCAGGTTCACCTAC")
primer_set_bx = PrimerSet(name="bx", fwd="GCCAGTAGTCATATGCTTGTCT", rev="GCCTGCTGCCTTCCTT")
primer_set_gf = PrimerSet(name="gf", fwd="", rev="")
primer_set_hd = PrimerSet(name="hd", fwd="GGACGATAAGACCCTATAAA", rev="ACGCTGTTATCCCTAAAGT")
primer_set_lm = PrimerSet(name="lm", fwd="", rev="")
primer_set_rv = PrimerSet(name="rv", fwd="TAGAACAGGCTCCTCTAG", rev="TTAGATACCCCACTATGC")
primer_set_um = PrimerSet(name="um", fwd="GGATTAGATACCCTGGTA", rev="CCGTCAATTCMTTTRAGTTT")
primer_set_wv = PrimerSet(name="wv", fwd="GACGAGAAGACCCTWTGGAGC", rev="CCRYGGTCGCCCCAAC")

database = ReferenceDatabase(dataset_ncbi_coi, primer_set_ci)
builder = DatabaseBuilder(database, working_dir="/data/pieter/workdir", environment=None, dry_run=False)
builder.build()
