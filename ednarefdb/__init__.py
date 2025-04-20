import os
import subprocess
import logging
from dataclasses import dataclass
from ednarefdb.lineage2taxtrain import lineage2taxtrain
from ednarefdb.addfulllineage import addfulllineage
import glob
from Bio import SeqIO
from ednarefdb.taxonomy import tsv_to_filled_taxonomy_tsv, fix_homonyms, generate_fasta_and_taxonomy


@dataclass
class PrimerSet:
    name: str
    fwd: str
    rev: str


@dataclass
class NucleotideDataset:
    name: str
    query: str = None
    path: str = None
    crabs_path: str = None


@dataclass
class ReferenceDatabase:
    dataset: NucleotideDataset
    primer_set: PrimerSet


class DatabaseBuilder:

    def __init__(self, working_dir: str, environment: str = "conda"):
        self.working_dir = working_dir
        self.environment = environment
        self.set_shell_config()

        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)

    def set_shell_config(self):
        shell = os.environ.get("SHELL")
        if shell:
            shell_name = os.path.basename(shell)
            if shell_name == "zsh":
                self.shell_config = "~/.zshrc"
                self.shell_executable = "/bin/zsh"
            elif shell_name == "bash":
                self.shell_config = "~/.bashrc"
                self.shell_executable = "/bin/bash"
            else:
                raise ValueError(f"Unsupported shell: {shell_name}")
        else:
            raise RuntimeError("Could not determine shell")

    def run_command_base(self, command):
        logging.info(command)
        subprocess.run(command, cwd=self.working_dir, shell=True, executable=self.shell_executable)

    # def run_command_conda(self, command):
    #     logging.info(command)
    #     subprocess.run(f"source {self.shell_config} && conda init && conda activate crabs && {command}", cwd=self.working_dir, shell=True, executable=self.shell_executable)

    def run_command_conda(self, command):
        # full_command = f"conda run -n crabs {command}"
        full_command = f"source {self.shell_config} && source $(conda info --base)/etc/profile.d/conda.sh && conda activate crabs && {command}"
        logging.info(full_command)
        subprocess.run(
            full_command,
            cwd=self.working_dir, 
            shell=True, 
            executable=self.shell_executable
        )

    def run_command_docker(self, command):
        logging.info(command)
        subprocess.run(f"docker run --rm -it -v $(pwd):/data --workdir='/data' quay.io/swordfish/crabs:0.1.4 {command}", cwd=self.working_dir, shell=True, executable=self.shell_executable)

    def run_command(self, command):
        if self.environment == "conda":
            self.run_command_conda(command)
        elif self.environment == "docker":
            self.run_command_docker(command)
        else:
            self.run_command_base(command)

    def cleanup_files(self, files: list[str]):
        for file in files:
                for f in glob.glob(os.path.join(self.working_dir, file)):
                    logging.info(f"Removing file: {f}")
                    try:
                        os.remove(f)
                    except FileNotFoundError:
                        pass

    def ncbi_download_taxonomy(self):

        logging.info(f"Downloading NCBI taxonomy to {self.working_dir}")
        self.run_command("crabs --download-taxonomy")

    def ncbi_download_nucleotide(self, dataset: NucleotideDataset):

        logging.info(f"Downloading NCBI nucleotide with query {dataset.query} to {dataset.name}")

        self.cleanup_files([
            f"{dataset.name}.fasta"
        ])

        self.run_command(f"""
            crabs --download-ncbi \
            --database nucleotide \
            --query '{dataset.query}' \
            --output {dataset.name}.fasta \
            --email helpdesk@obis.org \
            --batchsize 5000
        """)

    def import_nucleotide(self, dataset: NucleotideDataset):

        logging.info(f"Importing fasta for {dataset.name}")

        # TODO: set filenames at module level
        dataset_path = f"{dataset.name}.fasta" if dataset.path is None else dataset.path
        output_path = f"{dataset.name}.txt" if dataset.crabs_path is None else dataset.crabs_path

        self.run_command(f"""
            crabs --import \
            --import-format ncbi \
            --input {dataset_path} \
            --names names.dmp \
            --nodes nodes.dmp \
            --acc2tax nucl_gb.accession2taxid \
            --output {output_path} \
            --ranks 'domain;phylum;class;order;family;genus;species'
        """)

    def pcr(self, dataset: NucleotideDataset, primer_set: PrimerSet):

        logging.info(f"Performing in silico PCR for {dataset.name}_{primer_set.name}")

        self.cleanup_files([
            f"{dataset.name}_{primer_set.name}.txt"
        ])

        dataset_path = f"{dataset.name}.txt"

        self.run_command(f"""
            crabs --in-silico-pcr \
            --input {dataset_path} \
            --output {dataset.name}_{primer_set.name}.txt \
            --forward {primer_set.fwd} \
            --reverse {primer_set.rev}
        """)

    def pga(self, dataset: NucleotideDataset, primer_set: PrimerSet, percid: float = 0.8, coverage: float = 0.8):

        logging.info(f"Performing PGA for {dataset.name}_{primer_set.name}")

        self.cleanup_files([
            f"{dataset.name}_{primer_set.name}_pga.txt"
        ])

        dataset_path = f"{dataset.name}.txt"

        self.run_command(f"""
            crabs --pairwise-global-alignment \
            --input {dataset_path} \
            --output {dataset.name}_{primer_set.name}_pga.txt \
            --amplicons {dataset.name}_{primer_set.name}.txt \
            --forward {primer_set.fwd} \
            --reverse {primer_set.rev} \
            --percent-identity {percid} \
            --coverage {coverage}
        """)

    def dereplicate(self, dataset: NucleotideDataset, primer_set: PrimerSet):

        logging.info(f"Performing cleanup for {dataset.name}_{primer_set.name}")

        self.cleanup_files([
            f"{dataset.name}_{primer_set.name}_pga_derep.txt",
            # f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean.tsv",
            # f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.fasta",
            # f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_sintax.fasta",
            # f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.tsv",
            # f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled.tsv",
            # f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled_nona.tsv"
        ])

        self.run_command(f"""
            crabs --dereplicate \
            --input {dataset.name}_{primer_set.name}_pga.txt \
            --output {dataset.name}_{primer_set.name}_pga_derep.txt \
            --dereplication-method 'unique_species'
        """)

    def export(self, dataset: NucleotideDataset, primer_set: PrimerSet):

        logging.info(f"Performing export for {dataset.name}_{primer_set.name}")

        self.cleanup_files([
            f"{dataset.name}_{primer_set.name}_pga_derep_sintax.fasta",
            f"{dataset.name}_{primer_set.name}_pga_derep_rdp.fasta",
            f"{dataset.name}_{primer_set.name}_pga_derep_blast*"
        ])

        self.run_command(f"""
            crabs --export --export-format 'sintax' \
            --input {dataset.name}_{primer_set.name}_pga_derep.txt \
            --output {dataset.name}_{primer_set.name}_pga_derep_sintax.fasta
        """)

    def prepare_train(self, dataset: NucleotideDataset, primer_set: PrimerSet):

        logging.info(f"Preparing training files for {dataset.name}_{primer_set.name}")

        self.cleanup_files([
            # f"ready4train_{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.fasta",
            # f"ready4train_{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_names.fasta",
            # "Seq_IDs.tsv",
            # f"ready4train_{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_names2.fasta"
        ])

        input_file = os.path.join(self.working_dir, f"{dataset.name}_{primer_set.name}_pga_derep.txt")
        filled_file = os.path.join(self.working_dir, f"{dataset.name}_{primer_set.name}_pga_derep_filled.txt")
        fixed_file = os.path.join(self.working_dir, f"{dataset.name}_{primer_set.name}_pga_derep_filled_fixed.txt")
        rdp_tax_file = os.path.join(self.working_dir, f"{dataset.name}_{primer_set.name}_pga_derep_filled_fixed_rdp.txt")
        rdp_fasta_file = os.path.join(self.working_dir, f"{dataset.name}_{primer_set.name}_pga_derep_filled_fixed_rdp.fasta")
        rdp_taxonomy_file = os.path.join(self.working_dir, f"{dataset.name}_{primer_set.name}_pga_derep_filled_fixed_rdp_taxonomy.txt")

        tsv_to_filled_taxonomy_tsv(input_file, filled_file)
        fix_homonyms(filled_file, fixed_file)
        generate_fasta_and_taxonomy(fixed_file, rdp_fasta_file, rdp_tax_file)
        lineage2taxtrain(rdp_tax_file, rdp_taxonomy_file)

    def train(self, dataset: NucleotideDataset, primer_set: PrimerSet):

        logging.info(f"Training classifier for {dataset.name}_{primer_set.name}")

        rdp_fasta_file = os.path.join(self.working_dir, f"{dataset.name}_{primer_set.name}_pga_derep_filled_fixed_rdp.fasta")
        rdp_taxonomy_file = os.path.join(self.working_dir, f"{dataset.name}_{primer_set.name}_pga_derep_filled_fixed_rdp_taxonomy.txt")

        self.run_command_conda(f"""
            classifier -Xmx100g train -o training_files_{primer_set.name} \
            -s {rdp_fasta_file} \
            -t {rdp_taxonomy_file}
        """)

        # with open(os.path.join(self.working_dir, f"training_files_{primer_set.name}", "rRNAClassifier.properties"), "w") as xml_file:
        #     xml_file.write("bergeyTree=bergeyTrainingTree.xml" + "\n")
        #     xml_file.write("probabilityList=genus_wordConditionalProbList.txt" + "\n")
        #     xml_file.write("probabilityIndex=wordConditionalProbIndexArr.txt" + "\n")
        #     xml_file.write("wordPrior=logWordPrior.txt" + "\n")
        #     xml_file.write("classifierVersion=RDP Naive Bayesian rRNA Classifier Version ?")

    def cleanup(self, dataset: NucleotideDataset, primer_set: PrimerSet):

        self.cleanup_files([
            f"{dataset.name}_{primer_set.name}.fasta",
            f"{dataset.name}_{primer_set.name}_pga.fasta",
            f"{dataset.name}_{primer_set.name}_pga_taxa.tsv",
            f"{dataset.name}_{primer_set.name}_pga_missing_taxa.tsv",
            f"{dataset.name}_{primer_set.name}_pga_pacman.tsv",
            f"{dataset.name}_{primer_set.name}_pga_taxa_pacmanformat.tsv",
            f"{dataset.name}_{primer_set.name}_pga_taxa_derep.tsv",
            f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean.tsv",
            f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.fasta",
            f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.tsv",
            f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_names.fasta",
            f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_names2.fasta",
            f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled.tsv",
            f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled_nona.tsv",
            f"ready4train_{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.fasta",
            f"ready4train_{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.tsv",
            f"ready4train_{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_names.fasta",
            "Seq_IDs.tsv",
            f"ready4train_{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_names2.fasta"
        ])
