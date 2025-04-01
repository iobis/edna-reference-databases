from ednarefdb import PrimerSet, DatabaseBuilder, NucleotideDataset, ReferenceDatabase
import logging


logging.basicConfig(level=logging.INFO)


primer_set_mifish = PrimerSet(name="mifish", fwd="GTYGGTAAAWCTCGTGCCAGC", rev="CATAGTGGGGTATCTAATCCYAGTTTG")
dataset_nt = NucleotideDataset(
    name="ncbi_nt_filtered",
    path="/Volumes/acasis/nt_filtered/nt_filtered.fasta",
    crabs_path="/Volumes/acasis/nt_filtered/nt_filtered.crabs.txt"
)
databases = [
    ReferenceDatabase(dataset_nt, primer_set_mifish)
]


def main():

    builder = DatabaseBuilder(working_dir="./output", environment="conda")

    # builder.ncbi_download_taxonomy()
    # builder.ncbi_download_nucleotide(dataset)
    builder.import_nucleotide(dataset_nt)

    # for database in databases:
    #     print(database)
        # builder.pcr(dataset=database.dataset, primer_set=database.primer_set)
        # builder.pga(dataset=database.dataset, primer_set=database.primer_set, percid=0.8, coverage=0.8)
        # builder.assign_taxonomy(dataset=database.dataset, primer_set=database.primer_set)
        # builder.sequence_cleanup(dataset=database.dataset, primer_set=database.primer_set)
        # builder.train(dataset=database.dataset, primer_set=database.primer_set)
        # builder.cleanup(dataset=database.dataset, primer_set=database.primer_set)


if __name__ == "__main__":
    main()
