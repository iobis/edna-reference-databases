import os
import subprocess
import logging
from dataclasses import dataclass
from ednarefdb.lineage2taxtrain import lineage2taxtrain
from ednarefdb.addfulllineage import addfulllineage


@dataclass
class PrimerSet:
    name: str
    fwd: str
    rev: str


@dataclass
class NucleotideDataset:
    name: str
    query: str


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
            elif shell_name == "bash":
                self.shell_config = "~/.bashrc"
            else:
                raise ValueError(f"Unsupported shell: {shell_name}")
        else:
            raise RuntimeError("Could not determine shell")

    def run_command_conda(self, command):
        logging.info(command)
        subprocess.run(f"source {self.shell_config} && conda activate crabs && {command}", cwd=self.working_dir, shell=True, executable="/bin/bash")

    def run_command_docker(self, command):
        logging.info(command)
        subprocess.run(f"docker run --rm -it -v $(pwd):/data --workdir='/data' quay.io/swordfish/crabs:0.1.4 {command}", cwd=self.working_dir, shell=True, executable="/bin/bash")

    def run_command(self, command):
        if self.environment == "conda":
            self.run_command_conda(command)
        elif self.environment == "docker":
            self.run_command_docker(command)
        else:
            raise ValueError("Invalid environment")

    def cleanup_files(self, files: list[str]):
        for file in files:
            try:
                os.remove(os.path.join(self.working_dir, file))
            except FileNotFoundError:
                pass

    def ncbi_download_taxonomy(self):

        logging.info(f"Downloading NCBI taxonomy to {self.working_dir}")
        self.run_command("crabs db_download --source taxonomy")

    def ncbi_download_nucleotide(self, dataset: NucleotideDataset):

        logging.info(f"Downloading NCBI nucleotide with query {dataset.query} to {dataset.name}")

        self.cleanup_files([
            f"{dataset.name}.fasta"
        ])

        self.run_command(f"""
            crabs db_download \
            --source ncbi \
            --database nucleotide \
            --query '{dataset.query}' \
            --output {dataset.name}.fasta \
            --keep_original yes \
            --email helpdesk@obis.org \
            --batchsize 5000
        """)

    def pcr(self, dataset: str, primer_set: PrimerSet):

        logging.info(f"Performing in silico PCR for {dataset.name}_{primer_set.name}")

        self.cleanup_files([
            f"{dataset.name}_{primer_set.name}.fasta"
        ])

        self.run_command(f"""
            crabs insilico_pcr \
            --input {dataset.name}.fasta \
            --output {dataset.name}_{primer_set.name}.fasta \
            --fwd {primer_set.fwd} \
            --rev {primer_set.rev} \
            --error 4.5
        """)

    def pga(self, dataset: NucleotideDataset, primer_set: PrimerSet, percid: float = 0.8, coverage: float = 0.8):

        logging.info(f"Performing PGA for {dataset.name}_{primer_set.name}")

        self.cleanup_files([
            f"{dataset.name}_{primer_set.name}_pga.fasta"
        ])

        self.run_command(f"""
            crabs pga \
            --input {dataset.name}.fasta \
            --output {dataset.name}_{primer_set.name}_pga.fasta \
            --database {dataset.name}_{primer_set.name}.fasta \
            --fwd {primer_set.fwd} \
            --rev {primer_set.rev} \
            --speed medium \
            --percid {percid} \
            --coverage {coverage} \
            --filter_method relaxed
        """)

    def assign_taxonomy(self, dataset: NucleotideDataset, primer_set: PrimerSet):

        logging.info(f"Assigning taxonomy for {dataset.name}_{primer_set.name}")

        self.cleanup_files([
            f"{dataset.name}_{primer_set.name}_pga_taxa.tsv",
            f"{dataset.name}_{primer_set.name}_pga_missing_taxa.tsv",
            f"{dataset.name}_{primer_set.name}_pga_pacman.tsv",
            f"{dataset.name}_{primer_set.name}_pga_taxa_pacmanformat.tsv"
        ])

        self.run_command(f"""
            crabs assign_tax \
            --input {dataset.name}_{primer_set.name}_pga.fasta \
            --output {dataset.name}_{primer_set.name}_pga_taxa.tsv \
            --acc2tax nucl_gb.accession2taxid \
            --taxid nodes.dmp \
            --name names.dmp \
            --missing {dataset.name}_{primer_set.name}_pga_missing_taxa.tsv
        """)

        self.run_command(f"""
            awk -F'\t' '{{OFS="\t"; $10=$3";"$4";"$5";"$6";"$7";"$8";"$9; print}}' {dataset.name}_{primer_set.name}_pga_taxa.tsv > {dataset.name}_{primer_set.name}_pga_pacman.tsv
        """)

        self.run_command(f"""
            cut -f1,10 {dataset.name}_{primer_set.name}_pga_pacman.tsv > {dataset.name}_{primer_set.name}_pga_taxa_pacmanformat.tsv
        """)

    def sequence_cleanup(self, dataset: NucleotideDataset, primer_set: PrimerSet):

        logging.info(f"Performing cleanup for {dataset.name}_{primer_set.name}")

        self.cleanup_files([
            f"{dataset.name}_{primer_set.name}_pga_taxa_derep.tsv",
            f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean.tsv",
            f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.fasta",
            f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_sintax.fasta",
            f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.tsv",
            f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled.tsv",
            f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled_nona.tsv"
        ])

        self.run_command(f"""
            crabs dereplicate \
            --input {dataset.name}_{primer_set.name}_pga_taxa.tsv \
            --output {dataset.name}_{primer_set.name}_pga_taxa_derep.tsv \
            --method uniq_species
        """)

        self.run_command(f"""
            crabs seq_cleanup \
            --input {dataset.name}_{primer_set.name}_pga_taxa_derep.tsv \
            --maxns 0 \
            --output {dataset.name}_{primer_set.name}_pga_taxa_derep_clean.tsv \
            --minlen 100 --maxlen 5000 --nans 6 --enviro yes --species yes
        """)

        self.run_command(f"""
            awk '!/(\t|^)nan(\t|$)/' {dataset.name}_{primer_set.name}_pga_taxa_derep_clean.tsv > temp.tsv && mv temp.tsv {dataset.name}_{primer_set.name}_pga_taxa_derep_clean.tsv
        """)

        self.run_command(f"""
            crabs tax_format \
            --input {dataset.name}_{primer_set.name}_pga_taxa_derep_clean.tsv \
            --output {dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.fasta \
            --format rdp
        """)

        self.run_command(f"""
            crabs tax_format \
            --input {dataset.name}_{primer_set.name}_pga_taxa_derep_clean.tsv \
            --output {dataset.name}_{primer_set.name}_pga_taxa_derep_clean_sintax.fasta \
            --format sintax
        """)

        self.run_command(f"""
            grep "^>" {dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.fasta | \
            sed 's/^>//g' | \
            sed 's/;/\t/' \
            > {dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.tsv
        """)

        self.run_command(f"""
            awk 'BEGIN{{OFS=";"}} {{
            split($3, tax, ";");
            for (i = 1; i <= 7; i++) {{
                if (tax[i] == "") {{
                tax[i] = substr("kpcofgs", i, 1) "_" tax[i - 1];
                }}
            }}
            print $1, tax[1], tax[2], tax[3], tax[4], tax[5], tax[6], tax[7];
            }}' \
            {dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.tsv | \
            sed 's/;/\\t/g '> \
            {dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled.tsv
        """)

        self.run_command(f"""
            grep -v "g_f_o_c_p_k_" \
            {dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled.tsv > \
            {dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled_nona.tsv
        """)

        self.run_command(f"""
            wc -l {dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled_nona.tsv
        """)

        with open(os.path.join(self.working_dir, f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled_nona.tsv"), "r") as file:
            content = file.read()
        with open(os.path.join(self.working_dir, f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled_nona.tsv"), "w") as file:
            file.write("Seq-ID\tSuperkingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n")
            file.write(content)

    def train(self, dataset: NucleotideDataset, primer_set: PrimerSet):

        logging.info(f"Training classifier for {dataset.name}_{primer_set.name}")

        self.cleanup_files([
            f"ready4train_{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.fasta",
            f"ready4train_{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_names.fasta",
            "Seq_IDs.tsv",
            f"ready4train_{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_names2.fasta"
        ])

        lineage2taxtrain(
            os.path.join(self.working_dir, f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled_nona.tsv"),
            os.path.join(self.working_dir, f"ready4train_{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.tsv")
        )

        self.run_command_conda(f"""
            cut -f1 {dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.fasta > {dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_names.fasta
        """)

        self.run_command_conda(f"""
            cut -f1 {dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled_nona.tsv > Seq_IDs.tsv
        """)

        self.run_command_conda(f"""
            seqtk subseq {dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_names.fasta Seq_IDs.tsv > {dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_names2.fasta
        """)

        self.run_command_conda(f"""
            grep -c "^>" {dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_names2.fasta
        """)

        addfulllineage(
            os.path.join(self.working_dir, f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled_nona.tsv"),
            os.path.join(self.working_dir, f"{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp_names2.fasta"),
            os.path.join(self.working_dir, f"ready4train_{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.fasta")
        )

        self.run_command_conda(f"""
            classifier -Xmx16g train -o training_files_{primer_set.name} \
            -s ready4train_{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.fasta \
            -t ready4train_{dataset.name}_{primer_set.name}_pga_taxa_derep_clean_rdp.tsv
        """)

        with open(os.path.join(self.working_dir, f"training_files_{primer_set.name}", "rRNAClassifier.properties"), "w") as xml_file:
            xml_file.write("bergeyTree=bergeyTrainingTree.xml" + "\n")
            xml_file.write("probabilityList=genus_wordConditionalProbList.txt" + "\n")
            xml_file.write("probabilityIndex=wordConditionalProbIndexArr.txt" + "\n")
            xml_file.write("wordPrior=logWordPrior.txt" + "\n")
            xml_file.write("classifierVersion=RDP Naive Bayesian rRNA Classifier Version ?")

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
