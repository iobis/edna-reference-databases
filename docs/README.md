
# Reference databases for the PacMAN pipeline  

The PacMAN pipeline uses a Bowtie-2 - BLCA algorithm for the first level of taxonomic assignment with a reference database. This works best if the reference database is cut to the region of interest. Therefore to run the eDNA expeditions samples, a reference database was built for each target loci, using the [CRABS creating reference databases workflow](https://github.com/gjeunen/reference_database_creator). 

The reference database creator can download sequences directly from NCBI based on a word search, or alternatively can directly download the mitofish database. The workflow consists of 

1. Downloading the ncbi taxonomic database 
2. Downloading the target reference database based on a word search
3. Searching for the target primers in the downloaded sequences
4. Aligning the sequences that had primers with the remaining reference database (because many sequences will not have the primer region available)
5. Building the taxonomic reference file based on the sequences collected

This document shows the searches and analyses made in March 2023 to build the reference databases.

## Script

Execute `run.sh` to run the entire workflow.

## Set up 

The crabs workflow was downloaded based on the instructions in Manual Installation section.  

The conda environment file listing the dependencies downloaded to run crabs:

```
name: basta_py3
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - _libgcc_mutex=0.1=conda_forge
  - _openmp_mutex=4.5=2_gnu
  - bzip2=1.0.8=h7f98852_4
  - c-ares=1.18.1=h7f98852_0
  - ca-certificates=2022.12.7=ha878542_0
  - curl=7.88.1=hdc1c0ab_0
  - gettext=0.21.1=h27087fc_0
  - keyutils=1.6.1=h166bdaf_0
  - krb5=1.20.1=h81ceb04_0
  - krona=2.8.1=pl5321hdfd78af_1
  - ld_impl_linux-64=2.40=h41732ed_0
  - leveldb=1.23=h0b3b9e0_1
  - libcurl=7.88.1=hdc1c0ab_0
  - libedit=3.1.20191231=he28a2e2_2
  - libev=4.33=h516909a_1
  - libffi=3.4.2=h7f98852_5
  - libgcc-ng=12.2.0=h65d4601_19
  - libgomp=12.2.0=h65d4601_19
  - libidn2=2.3.4=h166bdaf_0
  - libnghttp2=1.51.0=hff17c54_0
  - libnsl=2.0.0=h7f98852_0
  - libsqlite=3.40.0=h753d276_0
  - libssh2=1.10.0=hf14f497_3
  - libstdcxx-ng=12.2.0=h46fd767_19
  - libunistring=0.9.10=h7f98852_0
  - libuuid=2.32.1=h7f98852_1000
  - libzlib=1.2.13=h166bdaf_4
  - muscle=3.8.31=0
  - ncurses=6.3=h27087fc_1
  - openssl=3.0.8=h0b41bf4_0
  - perl=5.32.1=2_h7f98852_perl5
  - pip=23.0.1=pyhd8ed1ab_0
  - plyvel=1.5.0=py311h4dd048b_1
  - python=3.11.0=he550d4f_1_cpython
  - python_abi=3.11=3_cp311
  - readline=8.1.2=h0f457ee_0
  - setuptools=67.4.0=pyhd8ed1ab_0
  - snappy=1.1.9=hbd366e4_2
  - tk=8.6.12=h27826a3_0
  - tzdata=2022g=h191b570_0
  - wheel=0.38.4=pyhd8ed1ab_0
  - xz=5.2.6=h166bdaf_0
  - zlib=1.2.13=h166bdaf_4
  - pip:
      - argparse==1.4.0
      - basta==1.4
      - biopython==1.81
      - contourpy==1.0.7
      - cycler==0.11.0
      - fonttools==4.38.0
      - kiwisolver==1.4.4
      - matplotlib==3.7.0
      - matplotlip==0.2
      - numpy==1.24.2
      - packaging==23.0
      - pandas==1.5.3
      - pillow==9.4.0
      - pyparsing==3.0.9
      - python-dateutil==2.8.2
      - pytz==2022.7.1
      - six==1.16.0
      - tqdm==4.64.1
      - wget==3.2
```

## Biomarkers

The different loci analysed in the eDNA expeditions project and the reference databases and searches done for these are:

|Loci |F-primer| R-primer | Reference database | Search text |
|-----|--------|----------|--| ---|
|Mifish-UE| GT**Y**GGTAAA**W**CTCGTGCCAGC      |CATAGTGGGGTATCTAATCC**Y**AGTTTG | NCBI |12S_ncbi_mifish_pga.fasta | 
|Mimammal-UEB  |  GG**RY**TGGT**H**AATTTCGTGCCAGC    | CATAGTG**R**GGTATCTAATC**Y**CAGTTTG | NCBI | 12S_ncbi_mimammals_pga.fasta |
| Teleo - 12S | ACACCGCCCGTCACTCT | CTTCCGGTACACTTACCATG | NCBI | 12S_ncbi_teleo_pga.fasta |
| Leray - COI | GGWACWGGWTGAACWGTWTAYCCYCC | TAIACYTCIGGRTGICCRAARAAYCA | MIDORI | no searches made MIDORI_UNIQ_SP_NUC_GB246_CO1_QIIME.fasta |
| Leray - COI | GGWACWGGWTGAACWGTWTAYCCYCC | TAIACYTCIGGRTGICCRAARAAYCA | NCBI | coi_leray_pga.fasta |
| Vert - 16S | AGACGAGAAGACCCYdTGGAGCTT | GATCCAACATCGAGGTCGTAA | NCBI | 16S_ncbi_vert_pga.fasta |

For Mifish and Mimammal the primers used for searching the reference databases and the sequences were a consensus primer from the different versions available. In bold are marked the non-specific basepairs due to building this consensus sequence. 

It is also possible to build the reference database using the mitofish database. In this case we used a search of NCBI.

## Sequence and taxonomy files

The setup steps consist of downloading the databases of choice and the taxonomic dump of ncbi.

### Taxonomy

```
crabs db_download --source taxonomy
```

### Sequences

Mifish and teleo (221,918 sequences collected)

```
./crabs db_download 
--source ncbi 
--database nucleotide 
--query '12S[All Fields] AND mitochondrial[All Fields] AND ("1"[SLEN] : "50000"[SLEN])'
--output 12S_mito_ncbi_1_50000.fasta 
--keep_original yes 
--email s@gmail.com 
--batchsize 5000
```

A query for mammals was made separately as a trial though the sequences could have been extracted also from the previous download. A total of 321,306 sequences were collected.

```
12S[All Fields] AND Mammalia[Organism] AND ("1"[SLEN] : "50000"[SLEN])  
```

And the same for 16S, the query resulted in about 445,298 sequences downloaded.

```
16S[All Fields] AND "Eukaryota"[Organism] AND ("1"[SLEN] : "50000"[SLEN]) AND mitochondrial[All Fields]
```

For COI, this resulted in 4,571,767 downloaded:

```
COI[All Fields] OR CO1[All Fields] OR cytochrome oxidase subunit I[All Fields] OR cytochrome oxidase subunit 1[All Fields] OR cytochrome c oxidase subunit I[All Fields] OR cytochrome c oxidase subunit 1[All Fields] OR COX1[All Fields] AND ("50"[SLEN] : "50000"[SLEN])
```

## Extract amplicon sequences

According to the github instructions the amplicon sequences were extracted with the following commands:

```
./crabs insilico_pcr 
--input 12S_mito_ncbi_1_50000.fasta 
--output 12S_ncbi_mifish.fasta 
--fwd GTYGGTAAAWCTCGTGCCAGC 
--rev CATAGTGGGGTATCTAATCCYAGTTTG 
--error 4.5
```

Pairwise global alignment to database is performed to find also target sequences were the primer region has been removed.

```
./crabs pga 
--input 12S_mito_ncbi_1_50000.fasta 
--output 12S_ncbi_mifish_pga.fasta 
--database 12S_ncbi_mifish.fasta 
--fwd GTYGGTAAAWCTCGTGCCAGC 
--rev CATAGTGGGGTATCTAATCCYAGTTTG 
--speed medium 
--percid 0.8 
--coverage 0.8 
--filter_method relaxed
```

And likewise for mimammal, teleo and 16S using the 16S database. The number of sequences remaining after the pcr are:

| Biomarker| ref db |downloaded | PCR | PGA | 
|---|---|---|---|---|
| Mifish-UE | ncbi |221,918 | 64,025 | 85,358 |
| Mimammal-UEB | ncbi | 321,306 | 140,778 | 174,006 |
| Teleo | ncbi | (same as mifish) | 62,158 | 89,796 |
| 16S | ncbi| 445,298 |  274,324 | 322,700|
| COI | ncbi | 4,571,767 | 383,232 | 4,552,087 |
| Mitofish-Mifish | mitofish | 796,239  | 16,786 | 31,441 |

Mitofish added here for comparison, but in the trial data analysis runs, the ncbi extracted data was used. 

## Assign taxonomy

Next the taxonomy was assigned using the downloaded ncbi taxonomy dump.

```
./crabs assign_tax 
--input 12S_ncbi_mifish_pga.fasta 
--output 12S_ncbi_mifish_pga_taxa.tsv 
--acc2tax /home/ubuntu/data/databases/ncbi/taxonomy20230224/nucl_gb.accession2taxid  
--taxid /home/ubuntu/data/databases/ncbi/taxonomy20230224/nodes.dmp 
--name /home/ubuntu/data/databases/ncbi/taxonomy20230224/names.dmp 
--missing 12S_ncbi_mifish_pga_missing_taxa.tsv
```

No other cleanup was done for the reference files, but is available in the crabs workflow. 

The PacMAN pipeline requires a simple text file with the sequence-id tab separated from the taxonomic information, which are separated by ",".

The taxonomy file was formatted for the PacMAN pipeline with the following commands:

```
awk -F'\t' '{OFS="\t", $10=$3","$4","$5","$6","$7","$8","$9, print}' 12S_ncbi_mifish_pga_taxa.tsv > 12S_ncbi_mifish_pga_pacman.tsv

cut -f1,10 12S_ncbi_mifish_pga_pacman.tsv > 12S_ncbi_mifish_pga_taxa_pacmanformat.tsv
```
