import os
import subprocess
import logging
from dataclasses import dataclass
from ednarefdb.lineage2taxtrain import lineage2taxtrain


@dataclass
class PrimerSet:
    name: str
    fwd: str
    rev: str


class ReferenceDatabase:

    def __init__(self, working_dir: str):
        self.working_dir = working_dir

        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)

    def run_command(self, command):
        logging.info(command)
        subprocess.run(f"source ~/.zshrc && conda activate crabs && {command}", cwd=self.working_dir, shell=True, executable="/bin/bash")

    def ncbi_download_taxonomy(self):

        logging.info(f"Downloading NCBI taxonomy to {self.working_dir}")
        self.run_command("crabs db_download --source taxonomy")

    def ncbi_download_nucleotide(self, query: str, name: str):

        logging.info(f"Downloading NCBI nucleotide with query {query} to {name}")

        self.run_command(f"""
            crabs db_download \
            --source ncbi \
            --database nucleotide \
            --query '{query}' \
            --output {name}.fasta \
            --keep_original yes \
            --email helpdesk@obis.org \
            --batchsize 5000
        """)

    def pcr(self, input: str, primer_set: PrimerSet):

        logging.info(f"Performing in silico PCR on {input} with primer {primer_set.name}")

        self.run_command(f"""
            crabs insilico_pcr \
            --input {input}.fasta \
            --output {input}_{primer_set.name}.fasta \
            --fwd {primer_set.fwd} \
            --rev {primer_set.rev} \
            --error 4.5
        """)

    def pga(self, input: str, primer_set: PrimerSet, percid: float = 0.8, coverage: float = 0.8):

        logging.info(f"Performing PGA on {input} with primer {primer_set.name}")

        self.run_command(f"""
            crabs pga \
            --input {input}.fasta \
            --output {input}_{primer_set.name}_pga.fasta \
            --database {input}_{primer_set.name}.fasta \
            --fwd {primer_set.fwd} \
            --rev {primer_set.rev} \
            --speed medium \
            --percid {percid} \
            --coverage {coverage} \
            --filter_method relaxed
        """)

    def assign_taxonomy(self, input: str, primer_set: PrimerSet):

        self.run_command(f"""
            crabs assign_tax \
            --input {input}_{primer_set.name}_pga.fasta \
            --output {input}_{primer_set.name}_pga_taxa.tsv \
            --acc2tax nucl_gb.accession2taxid \
            --taxid nodes.dmp \
            --name names.dmp \
            --missing {input}_{primer_set.name}_pga_missing_taxa.tsv
        """)

    def cleanup(self, input: str, primer_set: PrimerSet):

        self.run_command(f"""
            crabs dereplicate \
            --input {input}_{primer_set.name}_pga_taxa.tsv \
            --output {input}_{primer_set.name}_pga_taxa_derep.tsv \
            --method uniq_species
        """)

        self.run_command(f"""
            crabs seq_cleanup \
            --input {input}_{primer_set.name}_pga_taxa_derep.tsv \
            --maxns 0 \
            --output {input}_{primer_set.name}_pga_taxa_derep_clean.tsv \
            --minlen 100 --maxlen 5000 --nans 6 --enviro yes --species yes
        """)

        self.run_command(f"""
            crabs tax_format \
            --input {input}_{primer_set.name}_pga_taxa_derep_clean.tsv \
            --output {input}_{primer_set.name}_pga_taxa_derep_clean_rdp.fasta \
            --format rdp
        """)

        self.run_command(f"""
            crabs tax_format \
            --input {input}_{primer_set.name}_pga_taxa_derep_clean.tsv \
            --output {input}_{primer_set.name}_pga_taxa_derep_clean_sintax.fasta \
            --format sintax
        """)

        self.run_command(f"""
            grep "^>" {input}_{primer_set.name}_pga_taxa_derep_clean_rdp.fasta | \
            sed 's/^>//g' | \
            sed 's/;/\t/' \
            > {input}_{primer_set.name}_pga_taxa_derep_clean_rdp.tsv
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
            {input}_{primer_set.name}_pga_taxa_derep_clean_rdp.tsv | \
            sed 's/;/\\t/g '> \
            {input}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled.tsv
        """)

        self.run_command(f"""
            grep -v "g_f_o_c_p_k_" \
            {input}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled.tsv > \
            {input}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled_nona.tsv
        """)

        self.run_command(f"""
            wc -l {input}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled_nona.tsv
        """)

        with open(os.path.join(self.working_dir, f"{input}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled_nona.tsv"), "r") as file:
            content = file.read()
        with open(os.path.join(self.working_dir, f"{input}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled_nona.tsv"), "w") as file:
            file.write("Seq-ID\tSuperkingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n")
            file.write(content)

    def train(self, input: str, primer_set: PrimerSet):

        lineage2taxtrain(
            os.path.join(self.working_dir, f"{input}_{primer_set.name}_pga_taxa_derep_clean_rdp_filled_nona.tsv"),
            os.path.join(self.working_dir, f"ready4train_{input}_{primer_set.name}_pga_taxa_derep_clean_rdp.tsv")
        )
