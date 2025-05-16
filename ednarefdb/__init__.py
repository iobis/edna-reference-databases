import os
import subprocess
import logging
from dataclasses import dataclass
from ednarefdb.lineage2taxtrain import lineage2taxtrain
from ednarefdb.addfulllineage import addfulllineage
import glob
from Bio import SeqIO
from ednarefdb.taxonomy import tsv_to_filled_taxonomy_tsv, fix_homonyms, generate_fasta_and_taxonomy
import shutil


@dataclass
class PrimerSet:
    name: str
    fwd: str
    rev: str
    min_length: int = None
    max_length: int = None
    max_n: int = 1
    mismatch: int = 4
    pga_percid: float = 0.7


@dataclass
class NucleotideFile:
    path: str
    type: str


@dataclass
class NucleotideDataset:
    name: str
    query: str = None
    type: str = None
    # path: str = None
    # crabs_path: str = None
    files: list[NucleotideFile] = None


@dataclass
class ReferenceDatabase:
    dataset: NucleotideDataset
    primer_set: PrimerSet


class DatabaseBuilder:

    def __init__(self, database: ReferenceDatabase, working_dir: str, environment: str = "conda", dry_run: bool = False):
        self.database = database
        self.working_dir = working_dir
        self.environment = environment
        self.dry_run = dry_run
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
        subprocess.run(command, cwd=self.working_dir, shell=True, executable=self.shell_executable)

    # def run_command_conda(self, command):
    #     logging.info(command)
    #     subprocess.run(f"source {self.shell_config} && conda init && conda activate crabs && {command}", cwd=self.working_dir, shell=True, executable=self.shell_executable)

    def run_command_conda(self, command):
        full_command = f"source {self.shell_config} && source $(conda info --base)/etc/profile.d/conda.sh && conda activate crabs && {command}"
        logging.info(full_command)
        subprocess.run(
            full_command,
            cwd=self.working_dir, 
            shell=True, 
            executable=self.shell_executable
        )

    def run_command_docker(self, command):
        full_command = f"docker run --rm -it -v $(pwd):/data --workdir='/data' quay.io/swordfish/crabs:0.1.4 {command}"
        logging.info(full_command)
        subprocess.run(full_command, cwd=self.working_dir, shell=True, executable=self.shell_executable)

    def run_command(self, command, force_conda=False):
        logging.info(command)
        if not self.dry_run:
            if self.environment == "conda" or force_conda:
                self.run_command_conda(command)
            elif self.environment == "docker":
                self.run_command_docker(command)
            else:
                self.run_command_base(command)

    def cleanup_files(self, files: list[str]):
        for file in files:
            for f in glob.glob(os.path.join(self.working_dir, file), recursive=True):
                logging.info(f"Removing: {f}")
                try:
                    if os.path.isfile(f) or os.path.islink(f):
                        os.remove(f)
                    elif os.path.isdir(f):
                        shutil.rmtree(f)
                except FileNotFoundError:
                    pass
                except Exception as e:
                    logging.warning(f"Failed to remove {f}: {e}")

    def ncbi_download_taxonomy(self):

        logging.info(f"Downloading NCBI taxonomy to {self.working_dir}")
        self.run_command("crabs --download-taxonomy")

    def ncbi_download_nucleotide(self, dataset: NucleotideDataset):

        logging.info(f"Downloading NCBI nucleotide with query {dataset.query} to {dataset.name}")

        output_path = f"{dataset.name}.fasta"
        self.cleanup_files([
            output_path
        ])

        self.run_command(f"""
            crabs --download-ncbi \
            --database nucleotide \
            --query '{dataset.query}' \
            --output {output_path} \
            --email helpdesk@obis.org \
            --batchsize 5000
        """)

    def import_nucleotide(self):

        logging.info(f"Importing fasta for {self.database.dataset.name}")

        output_path = f"{self.database.dataset.name}.txt"
        if not self.dry_run: self.cleanup_files([ output_path ])

        if self.database.dataset.files is not None:
            temp_files = []

            for i, file in enumerate(self.database.dataset.files):
                temp_output_path = f"{self.database.dataset.name}_{i}.txt"
                temp_files.append(temp_output_path)
                self.run_command(f"""
                    crabs --import \
                    --import-format {file.type} \
                    --input {file.path} \
                    --names names.dmp \
                    --nodes nodes.dmp \
                    --acc2tax nucl_gb.accession2taxid \
                    --output {temp_output_path} \
                    --ranks 'domain;phylum;class;order;family;genus;species'
                """)

            if len(temp_files) > 1:
                temp_files_concat = ";".join(temp_files)
                self.run_command(f"""
                    crabs --merge \
                    --input '{temp_files_concat}' \
                    --uniq \
                    --output {output_path}
                """)
            else:
                if not self.dry_run:
                    os.rename(os.path.join(self.working_dir, temp_files[0]), os.path.join(self.working_dir, output_path))

        else:
            input_path = f"{self.database.dataset.name}.fasta"
            self.run_command(f"""
                crabs --import \
                --input {input_path} \
                --import-format {self.database.dataset.type} \
                --names names.dmp \
                --nodes nodes.dmp \
                --acc2tax nucl_gb.accession2taxid \
                --output {output_path} \
                --ranks 'domain,phylum,class,order,family,genus,species'
            """)

    def pcr(self):

        logging.info(f"Performing in silico PCR for {self.database.dataset.name}_{self.database.primer_set.name}")

        dataset_path = f"{self.database.dataset.name}.txt"
        prefilter_path = f"{self.database.dataset.name}_{self.database.primer_set.name}_prefilter.txt"
        output_path = f"{self.database.dataset.name}_{self.database.primer_set.name}.txt"
        if not self.dry_run:
            self.cleanup_files([
                prefilter_path,
                output_path
            ])

        self.run_command(f"""
            crabs --in-silico-pcr \
            --input {dataset_path} \
            --threads 4 \
            --mismatch {self.database.primer_set.mismatch} \
            --output {prefilter_path} \
            --forward {self.database.primer_set.fwd} \
            --reverse {self.database.primer_set.rev}
        """)

        self.run_command(f"""
            crabs --filter \
            --input {prefilter_path} \
            --output {output_path} \
            --minimum-length {self.database.primer_set.min_length} \
            --maximum-length {self.database.primer_set.max_length} \
            --maximum-n {self.database.primer_set.max_n}
        """)


    def pga(self):

        logging.info(f"Performing PGA for {self.database.dataset.name}_{self.database.primer_set.name}")

        output_path = f"{self.database.dataset.name}_{self.database.primer_set.name}_pga.txt"
        if not self.dry_run: self.cleanup_files([ output_path ])

        dataset_path = f"{self.database.dataset.name}.txt"
        amplicon_path = f"{self.database.dataset.name}_{self.database.primer_set.name}.txt"

        self.run_command(f"""
            crabs --pairwise-global-alignment \
            --input {dataset_path} \
            --output {output_path} \
            --amplicons {amplicon_path} \
            --forward {self.database.primer_set.fwd} \
            --reverse {self.database.primer_set.rev} \
            --percent-identity {self.database.primer_set.pga_percid} \
            --coverage 1
        """)

    def dereplicate(self):

        logging.info(f"Performing cleanup for {self.database.dataset.name}_{self.database.primer_set.name}")

        input_path = f"{self.database.dataset.name}_{self.database.primer_set.name}_pga.txt"
        output_path = f"{self.database.dataset.name}_{self.database.primer_set.name}_pga_derep.txt"
        if not self.dry_run: self.cleanup_files([ output_path ])

        self.run_command(f"""
            crabs --dereplicate \
            --input {input_path} \
            --output {output_path} \
            --dereplication-method 'unique_species'
        """)

    def filter(self):

        logging.info(f"Performing cleanup for {self.database.dataset.name}_{self.database.primer_set.name}")

        input_path = "{self.database.dataset.name}_{self.database.primer_set.name}_pga_derep.txt"
        output_path = "{self.database.dataset.name}_{self.database.primer_set.name}_pga_derep_filtered.txt"
        if not self.dry_run: self.cleanup_files([ output_path ])

        self.run_command(f"""
            crabs --filter \
            --input {input_path} \
            --output {output_path} \
            --minimum-length {self.database.primer_set.min_length} \
            --maximum-length {self.database.primer_set.max_length} \
            --maximum-n {self.database.primer_set.max_n}
        """)

    def export(self):

        logging.info(f"Performing export for {self.database.dataset.name}_{self.database.primer_set.name}")

        input_path = "{self.database.dataset.name}_{self.database.primer_set.name}_pga_derep_filtered.txt"
        output_path = "{self.database.dataset.name}_{self.database.primer_set.name}_pga_derep_filtered_sintax.fasta"
        if not self.dry_run: self.cleanup_files([ output_path ])

        self.run_command(f"""
            crabs --export --export-format 'sintax' \
            --input {input_path} \
            --output {output_path}
        """)

    def prepare_train(self):

        logging.info(f"Preparing training files for {self.database.dataset.name}_{self.database.primer_set.name}")

        input_file = os.path.join(self.working_dir, f"{self.database.dataset.name}_{self.database.primer_set.name}_pga_derep_filtered.txt")
        filled_file = os.path.join(self.working_dir, f"{self.database.dataset.name}_{self.database.primer_set.name}_pga_derep_filtered_filled.txt")
        fixed_file = os.path.join(self.working_dir, f"{self.database.dataset.name}_{self.database.primer_set.name}_pga_derep_filtered_filled_fixed.txt")
        rdp_tax_file = os.path.join(self.working_dir, f"{self.database.dataset.name}_{self.database.primer_set.name}_pga_derep_filtered_filled_fixed_rdp.txt")
        rdp_fasta_file = os.path.join(self.working_dir, f"{self.database.dataset.name}_{self.database.primer_set.name}_pga_derep_filtered_filled_fixed_rdp.fasta")
        rdp_taxonomy_file = os.path.join(self.working_dir, f"{self.database.dataset.name}_{self.database.primer_set.name}_pga_derep_filtered_filled_fixed_rdp_taxonomy.txt")

        if not self.dry_run:
            self.cleanup_files([
                filled_file,
                fixed_file,
                rdp_tax_file,
                rdp_fasta_file,
                rdp_taxonomy_file
            ])

            tsv_to_filled_taxonomy_tsv(input_file, filled_file)
            fix_homonyms(filled_file, fixed_file)
            generate_fasta_and_taxonomy(fixed_file, rdp_fasta_file, rdp_tax_file)
            lineage2taxtrain(rdp_tax_file, rdp_taxonomy_file)

    def train(self):

        logging.info(f"Training classifier for {self.database.dataset.name}_{self.database.primer_set.name}")

        rdp_fasta_file = os.path.join(self.working_dir, f"{self.database.dataset.name}_{self.database.primer_set.name}_pga_derep_filtered_filled_fixed_rdp.fasta")
        rdp_taxonomy_file = os.path.join(self.working_dir, f"{self.database.dataset.name}_{self.database.primer_set.name}_pga_derep_filtered_filled_fixed_rdp_taxonomy.txt")
        output_path = f"training_files_{self.database.primer_set.name}"
        if not self.dry_run: self.cleanup_files([ output_path ])

        self.run_command(f"""
            classifier -Xmx100g train -o {output_path} \
            -s {rdp_fasta_file} \
            -t {rdp_taxonomy_file}
        """, force_conda=True)

        # with open(os.path.join(self.working_dir, f"training_files_{primer_set.name}", "rRNAClassifier.properties"), "w") as xml_file:
        #     xml_file.write("bergeyTree=bergeyTrainingTree.xml" + "\n")
        #     xml_file.write("probabilityList=genus_wordConditionalProbList.txt" + "\n")
        #     xml_file.write("probabilityIndex=wordConditionalProbIndexArr.txt" + "\n")
        #     xml_file.write("wordPrior=logWordPrior.txt" + "\n")
        #     xml_file.write("classifierVersion=RDP Naive Bayesian rRNA Classifier Version ?")

    def build(self):
        self.import_nucleotide()
        self.pcr()
        self.pga()
        self.dereplicate()
        self.filter()
        self.export()
        self.prepare_train()
        self.train()
