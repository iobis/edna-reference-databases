from ednarefdb import ReferenceDatabase, PrimerSet
import logging


logging.basicConfig(level=logging.INFO)


primer_set_12s = PrimerSet(
    name="mifish",
    fwd="GTYGGTAAAWCTCGTGCCAGC",
    rev="CATAGTGGGGTATCTAATCCYAGTTTG"
)


def main():

    ref_db = ReferenceDatabase(working_dir="./output")

    # ref_db.ncbi_download_taxonomy()

    # ref_db.ncbi_download_nucleotide(
    #     query='12S[All Fields] AND mitochondrial[All Fields] AND ("1"[SLEN] : "50000"[SLEN])',
    #     name="12s_mito_1_50000"
    # )

    # ref_db.pcr(
    #     input="12s_mito_1_50000",
    #     primer_set=primer_set_12s
    # )

    # ref_db.pga(
    #     input="12s_mito_1_50000",
    #     primer_set=primer_set_12s,
    #     percid=0.8,
    #     coverage=0.8
    # )

    # ref_db.assign_taxonomy(
    #     input="12s_mito_1_50000",
    #     primer_set=primer_set_12s
    # )

    # ref_db.cleanup(
    #     input="12s_mito_1_50000",
    #     primer_set=primer_set_12s
    # )

    ref_db.train(
        input="12s_mito_1_50000",
        primer_set=primer_set_12s
    )


if __name__ == "__main__":
    main()
