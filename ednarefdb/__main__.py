from ednarefdb import PrimerSet, DatabaseBuilder, NucleotideDataset, ReferenceDatabase
import logging


logging.basicConfig(level=logging.INFO)


primer_set_mifish = PrimerSet(name="mifish", fwd="GTYGGTAAAWCTCGTGCCAGC", rev="CATAGTGGGGTATCTAATCCYAGTTTG")
primer_set_mimammal = PrimerSet(name="mimammal", fwd="GGRYTGGTHAATTTCGTGCCAGC", rev="CATAGTGRGGTATCTAATCYCAGTTTG")
primer_set_teleo = PrimerSet(name="teleo", fwd="ACACCGCCCGTCACTCT", rev="CTTCCGGTACACTTACCATG")
primer_set_16s = PrimerSet(name="16s", fwd="AGACGAGAAGACCCYdTGGAGCTT", rev="GATCCAACATCGAGGTCGTAA")
primer_set_coi = PrimerSet(name="coi", fwd="GGWACWGGWTGAACWGTWTAYCCYCC", rev="TAIACYTCIGGRTGICCRAARAAYCA")

dataset_12s_ribo = NucleotideDataset(name="12s_ncbi_ribo_1_50000", query='12S[All Fields] AND ribosomal[All Fields] AND ("1"[SLEN] : "50000"[SLEN])')
dataset_16s_ribo = NucleotideDataset(name="16s_ncbi_ribo_1_50000", query='16S[All Fields] AND ribosomal[All Fields] AND ("1"[SLEN] : "50000"[SLEN])')
dataset_12s16s_ribo = NucleotideDataset(name="12s16s_ncbi_ribo_1_50000", query='(12S[All Fields] OR 16S[All Fields]) AND ribosomal[All Fields] AND ("1"[SLEN] : "50000"[SLEN])')
dataset_coi = NucleotideDataset(name="coi_ncbi_1_50000", query='COI[All Fields] OR CO1[All Fields] OR cytochrome oxidase subunit I[All Fields] OR cytochrome oxidase subunit 1[All Fields] OR cytochrome c oxidase subunit I[All Fields] OR cytochrome c oxidase subunit 1[All Fields] OR COX1[All Fields] AND ("50"[SLEN] : "50000"[SLEN])')

datasets = [
    # dataset_12s_ribo,
    # dataset_16s_ribo,
    # dataset_12s16s_ribo,
    dataset_coi
]

databases = [
    # ReferenceDatabase(dataset_12s16s_ribo, primer_set_mifish),
    # ReferenceDatabase(dataset_12s16s_ribo, primer_set_mimammal),
    # ReferenceDatabase(dataset_12s16s_ribo, primer_set_teleo),
    # ReferenceDatabase(dataset_16s_ribo, primer_set_16s),
    ReferenceDatabase(dataset_coi, primer_set_coi)
]


def main():

    ref_db = DatabaseBuilder(working_dir="./output", environment="conda")

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


if __name__ == "__main__":
    main()
