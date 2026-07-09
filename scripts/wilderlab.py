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


dataset_ncbi_coi = NucleotideDataset(name="ncbi_coi", files=[NucleotideFile(path="/Volumes/acasis/ncbi/coi_sequences.fasta", type="ncbi")])
dataset_ncbi_12s = NucleotideDataset(name="ncbi_12s", files=[NucleotideFile(path="/Volumes/acasis/ncbi/12s_sequences.fasta", type="ncbi")])
dataset_ncbi_18s = NucleotideDataset(name="ncbi_18s", files=[NucleotideFile(path="/Volumes/acasis/ncbi/18s_sequences.fasta", type="ncbi")])
dataset_ncbi_28s = NucleotideDataset(name="ncbi_28s", files=[NucleotideFile(path="/Volumes/acasis/ncbi/28s_sequences.fasta", type="ncbi")])
dataset_ncbi_its = NucleotideDataset(name="ncbi_its", files=[NucleotideFile(path="/Volumes/acasis/ncbi/its_sequences.fasta", type="ncbi")])
dataset_ncbi_16s_eukaryotic = NucleotideDataset(name="ncbi_16s_eukaryotic", files=[NucleotideFile(path="/Volumes/acasis/ncbi/16s_eukaryotic_sequences.fasta", type="ncbi")])
dataset_ncbi_16s_prokaryotic = NucleotideDataset(name="ncbi_16s_prokaryotic", files=[NucleotideFile(path="/Volumes/acasis/ncbi/16s_prokaryotic_sequences.fasta", type="ncbi")])
dataset_ncbi_16s_prokaryotic_slim = NucleotideDataset(name="ncbi_16s_prokaryotic_slim", files=[NucleotideFile(path="/Volumes/acasis/ncbi/16s_prokaryotic_slim_sequences.fasta", type="ncbi")])

primer_set_ci = primer_set_from_assay("CI", mismatch=0, pga_percid=0.7)
primer_set_lx = primer_set_from_assay("LX", mismatch=0, pga_percid=0.7)
primer_set_rv = primer_set_from_assay("RV", mismatch=0, pga_percid=0.7)
primer_set_lm = primer_set_from_assay("LM", mismatch=0, pga_percid=0.7)
primer_set_bu = primer_set_from_assay("BU", mismatch=0, pga_percid=0.7)
primer_set_be = primer_set_from_assay("BE", mismatch=0, pga_percid=0.7)
primer_set_bx = primer_set_from_assay("BX", mismatch=0, pga_percid=0.7)
primer_set_gf = primer_set_from_assay("GF", mismatch=0, pga_percid=0.7)
primer_set_gf.max_length = 193
primer_set_gf_relaxed = primer_set_from_assay("GF", mismatch=4, pga_percid=0.5)
primer_set_gf_relaxed.max_length = 193
primer_set_gf_relaxed_long = primer_set_from_assay("GF", mismatch=4, pga_percid=0.5)
primer_set_gf_relaxed.max_length = 300
primer_set_gd = primer_set_from_assay("GD", mismatch=0, pga_percid=0.7)
primer_set_wv = primer_set_from_assay("WV", mismatch=0, pga_percid=0.7)
primer_set_hd = primer_set_from_assay("HD", mismatch=0, pga_percid=0.7)
primer_set_wg = primer_set_from_assay("WG", mismatch=0, pga_percid=0.7)
primer_set_um = primer_set_from_assay("UM", mismatch=0, pga_percid=0.7)
primer_set_mc = primer_set_from_assay("MC", mismatch=0, pga_percid=0.7)
primer_set_ma = primer_set_from_assay("MA", mismatch=0, pga_percid=0.7)

database_ci = ReferenceDatabase(dataset_ncbi_coi, primer_set_ci)
database_lx = ReferenceDatabase(dataset_ncbi_12s, primer_set_lx)
database_rv = ReferenceDatabase(dataset_ncbi_12s, primer_set_rv)
database_lm = ReferenceDatabase(dataset_ncbi_12s, primer_set_lm)
database_bu = ReferenceDatabase(dataset_ncbi_18s, primer_set_bu)
database_be = ReferenceDatabase(dataset_ncbi_18s, primer_set_be)
database_bx = ReferenceDatabase(dataset_ncbi_18s, primer_set_bx)
database_gf = ReferenceDatabase(dataset_ncbi_its, primer_set_gf)
database_gf_relaxed = ReferenceDatabase(dataset_ncbi_its, primer_set_gf_relaxed, suffix="relaxed")
database_gf_relaxed_long = ReferenceDatabase(dataset_ncbi_its, primer_set_gf_relaxed_long, suffix="relaxed_long")
database_gd = ReferenceDatabase(dataset_ncbi_its, primer_set_gd)
database_wv = ReferenceDatabase(dataset_ncbi_16s_eukaryotic, primer_set_wv)
database_hd = ReferenceDatabase(dataset_ncbi_16s_eukaryotic, primer_set_hd)
database_wg = ReferenceDatabase(dataset_ncbi_16s_eukaryotic, primer_set_wg)
database_um = ReferenceDatabase(dataset_ncbi_16s_prokaryotic_slim, primer_set_um)
database_mc = ReferenceDatabase(dataset_ncbi_28s, primer_set_mc)
database_ma = ReferenceDatabase(dataset_ncbi_28s, primer_set_ma)

# builder.ncbi_download_taxonomy()

for database in [
    # database_ci,
    # database_rv,
    # database_lm,
    # database_bu,
    # database_be,
    # database_bx,
    # database_lx,
    # database_gf,
    # database_gf_relaxed,
    database_gf_relaxed_long,
    # database_gd,
    # database_wv,
    # database_hd,
    # database_wg,
    # database_mc,
    # database_ma,
    # database_um,
]:
    builder = DatabaseBuilder(database, working_dir="./workdir", environment=None, dry_run=False)
    builder.build()
    builder.compress_sintax()
