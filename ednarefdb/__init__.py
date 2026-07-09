import os
import subprocess
import logging
import gzip
import shutil
from datetime import date
from dataclasses import dataclass
from ednarefdb.ncbi_fasta import NcbiFastaBuilder, NcbiFastaBuildResult
from ednarefdb.report import BuildReportCollector
import glob
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
    suffix: str | None = None

    @property
    def build_stem(self) -> str:
        stem = f"{self.dataset.name}_{self.primer_set.name}"
        if self.suffix:
            stem = f"{stem}_{self.suffix}"
        return stem


class DatabaseBuilder:

    def __init__(self, database: ReferenceDatabase, working_dir: str, environment: str = "conda", dry_run: bool = False):
        self.database = database
        self.working_dir = working_dir
        self.environment = environment
        self.dry_run = dry_run

        self.set_shell_config()
        self.crabs_env = self._resolve_crabs_env()

        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)

        self.report_collector = BuildReportCollector(
            working_dir=self.working_dir,
            dry_run=self.dry_run,
            environment=self.environment,
        )

    def path(self, key: str) -> str:
        return os.path.join(self.working_dir, self.files[key])

    def set_file_paths(self):
        name = self.database.dataset.name
        stem = self.database.build_stem
        self.files = {
            "dataset_crabs_file": f"{name}.txt",
            "pcr_prefilter_file": f"{stem}_prefilter.txt",
            "pcr_file": f"{stem}.txt",
            "pga_file": f"{stem}_pga.txt",
            "dereplicate_file": f"{stem}_pga_derep.txt",
            "filtered_file": f"{stem}_pga_derep_filtered.txt",
            "sintax_file": f"{stem}_pga_derep_filtered_sintax.fasta",
        }

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

    def _crabs_available(self, env: dict[str, str]) -> bool:
        result = subprocess.run(
            "crabs --help",
            shell=True,
            executable=self.shell_executable,
            env=env,
            capture_output=True,
        )
        return result.returncode == 0

    def _resolve_crabs_env(self) -> dict[str, str]:
        env = os.environ.copy()
        if self._crabs_available(env):
            return env

        result = subprocess.run(
            ["pyenv", "whence", "crabs"],
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            pyenv_version = result.stdout.strip().splitlines()[0]
            env["PYENV_VERSION"] = pyenv_version
            if self._crabs_available(env):
                logging.info(f"Using crabs from pyenv {pyenv_version}")
                return env

        raise RuntimeError(
            "crabs not found on PATH. Install crabs, set environment='conda', "
            "or use a pyenv version that provides it."
        )

    def run_command_base(self, command):
        result = subprocess.run(
            command,
            cwd=self.working_dir,
            shell=True,
            executable=self.shell_executable,
            env=self.crabs_env,
        )
        if result.returncode != 0:
            raise RuntimeError(f"Command failed with exit code {result.returncode}")

    # def run_command_conda(self, command):
    #     logging.info(command)
    #     subprocess.run(f"source {self.shell_config} && conda init && conda activate crabs && {command}", cwd=self.working_dir, shell=True, executable=self.shell_executable)

    def run_command_conda(self, command):
        full_command = f"source {self.shell_config} && source $(conda info --base)/etc/profile.d/conda.sh && conda activate crabs && {command}"
        logging.info(full_command)
        result = subprocess.run(
            full_command,
            cwd=self.working_dir,
            shell=True,
            executable=self.shell_executable
        )
        if result.returncode != 0:
            raise RuntimeError(f"Command failed with exit code {result.returncode}")

    def run_command_docker(self, command):
        full_command = f"docker run --rm -it -v $(pwd):/data --workdir='/data' quay.io/swordfish/crabs:0.1.4 {command}"
        logging.info(full_command)
        result = subprocess.run(
            full_command, cwd=self.working_dir, shell=True, executable=self.shell_executable
        )
        if result.returncode != 0:
            raise RuntimeError(f"Command failed with exit code {result.returncode}")

    def run_command(self, command, force_conda=False):
        logging.info(command)
        self.report_collector.record_command(command)
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

        with self.report_collector.step("Import"):
            if not self.dry_run: self.cleanup_files([ self.files["dataset_crabs_file"] ])

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
                        --output {self.files["dataset_crabs_file"]}
                    """)
                else:
                    if not self.dry_run:
                        os.rename(
                            os.path.join(self.working_dir, temp_files[0]),
                            self.path("dataset_crabs_file"),
                        )

            else:
                input_path = f"{self.database.dataset.name}.fasta"
                self.run_command(f"""
                    crabs --import \
                    --input {input_path} \
                    --import-format {self.database.dataset.type} \
                    --names names.dmp \
                    --nodes nodes.dmp \
                    --acc2tax nucl_gb.accession2taxid \
                    --output {self.files["dataset_crabs_file"]} \
                    --ranks 'domain,phylum,class,order,family,genus,species'
                """)

    def pcr(self):

        logging.info(f"Performing in silico PCR for {self.database.build_stem}")

        with self.report_collector.step("In-silico PCR"):
            if not self.dry_run:
                self.cleanup_files([
                    self.files["pcr_prefilter_file"],
                    self.files["pcr_file"]
                ])

            self.run_command(f"""
                crabs --in-silico-pcr \
                --input {self.files['dataset_crabs_file']} \
                --threads 4 \
                --mismatch {self.database.primer_set.mismatch} \
                --output {self.files['pcr_prefilter_file']} \
                --forward {self.database.primer_set.fwd} \
                --reverse {self.database.primer_set.rev}
            """)

            self.run_command(f"""
                crabs --filter \
                --input {self.files['pcr_prefilter_file']} \
                --output {self.files['pcr_file']} \
                --minimum-length {self.database.primer_set.min_length} \
                --maximum-length {self.database.primer_set.max_length} \
                --maximum-n {self.database.primer_set.max_n}
            """)


    def pga(self):

        logging.info(f"Performing PGA for {self.database.build_stem}")

        with self.report_collector.step("PGA"):
            if not self.dry_run: self.cleanup_files([ self.files['pga_file'] ])

            self.run_command(f"""
                crabs --pairwise-global-alignment \
                --input {self.files['dataset_crabs_file']} \
                --output {self.files['pga_file']} \
                --amplicons {self.files['pcr_file']} \
                --forward {self.database.primer_set.fwd} \
                --reverse {self.database.primer_set.rev} \
                --percent-identity {self.database.primer_set.pga_percid} \
                --coverage 1
            """)

    def dereplicate(self):

        logging.info(f"Performing cleanup for {self.database.build_stem}")

        with self.report_collector.step("Dereplicate"):
            if not self.dry_run: self.cleanup_files([ self.files['dereplicate_file'] ])

            self.run_command(f"""
                crabs --dereplicate \
                --input {self.files['pga_file']} \
                --output {self.files['dereplicate_file']} \
                --dereplication-method 'unique_species'
            """)

    def filter(self):

        logging.info(f"Performing cleanup for {self.database.build_stem}")

        with self.report_collector.step("Filter"):
            if not self.dry_run: self.cleanup_files([ self.files["filtered_file"] ])

            self.run_command(f"""
                crabs --filter \
                --input {self.files['dereplicate_file']} \
                --output {self.files['filtered_file']} \
                --minimum-length {self.database.primer_set.min_length} \
                --maximum-length {self.database.primer_set.max_length} \
                --maximum-n {self.database.primer_set.max_n}
            """)

    def export(self):

        logging.info(f"Performing export for {self.database.build_stem}")

        with self.report_collector.step("SINTAX export"):
            if not self.dry_run: self.cleanup_files([ self.files["sintax_file"] ])

            self.run_command(f"""
                crabs --export --export-format 'sintax' \
                --input {self.files['filtered_file']} \
                --output {self.files['sintax_file']}
            """)

    def count_sequences(self, filename):
        result = subprocess.run(
            f'grep ">" "{os.path.join(self.working_dir, filename)}" | wc -l',
            shell=True,
            capture_output=True,
            text=True
        )
        return int(result.stdout.strip().split()[0])

    def count_lines(self, filename):
        result = subprocess.run(
            f'wc -l {os.path.join(self.working_dir, filename)}',
            shell=True,
            capture_output=True,
            text=True
        )
        return int(result.stdout.strip().split()[0])

    def summarize_files(self):
        for file in self.files.values():
            full_path = os.path.join(self.working_dir, file)
            if os.path.exists(full_path) and file.endswith(".fasta"):
                logging.info(f"Number of sequences in {full_path}: {self.count_sequences(file)}")
            elif os.path.exists(full_path) and file.endswith(".txt"):
                logging.info(f"Number of lines in {full_path}: {self.count_lines(file)}")

    def generate_report(self) -> str:
        self.report_collector.set_database(self.database)
        report_path = self.report_collector.write_report(self.files)
        logging.info(f"Build report written to {report_path}")
        return report_path

    def compress_sintax(self, datestamp: str | None = None) -> str:
        source = self.path("sintax_file")
        if not os.path.exists(source):
            raise FileNotFoundError(f"SINTAX output not found: {source}")

        stamp = datestamp or date.today().strftime("%Y%m%d")
        stem = self.files["sintax_file"].removesuffix(".fasta")
        output_name = f"{stem}_{stamp}.fasta.gz"
        output_path = os.path.join(self.working_dir, output_name)

        with open(source, "rb") as handle_in:
            with gzip.open(output_path, "wb") as handle_out:
                shutil.copyfileobj(handle_in, handle_out)

        logging.info(f"Compressed SINTAX output to {output_path}")
        return output_path

    def build(self):
        self.set_file_paths()
        self.report_collector.set_database(self.database)
        self.report_collector.begin_build()
        self.import_nucleotide()
        self.pcr()
        self.pga()
        self.dereplicate()
        self.filter()
        self.export()
        self.report_collector.finish_build()
        self.summarize_files()
        self.generate_report()
