# Training an rdp reference database

This document shows the workflow of training a reference database for the rdp classifier following the tutorial of [John Quensen](https://john-quensen.com/tutorials/training-the-rdp-classifier/). 

The databases are used for the classification of reads in the [eDNA expeditions project](https://www.unesco.org/en/edna-expeditions), and the assignment method will be integrated into the [PacMAN pipeline](https://github.com/iobis/PacMAN-pipeline). The species lists will be reviewed by local site managers on the [sample dashboard](https://samples.ednaexpeditions.org/).

The eDNA expeditions project uses five different biomarkers for the analysis of biodiversity from the samples. Reference databases for each of these different markers were build from NCBI using the [CRABS](https://github.com/gjeunen/reference_database_creator) workflow. 

As the first example, we train a reference database for the COI marker used in the project. The results are compared to a comprehensive trained database available from [T. M. Porter](https://github.com/terrimporter/CO1Classifier). They also provide a 12S database that is already trained. However eDNA expeditions includes also a 16S marker for fish, so training our own databases is necessary. 

## Download and build the database


### Workflow for the previous version of PacMAN

The PacMAN pipeline uses a Bowtie-2 - BLCA algorithm for the first level of taxonomic assignment with a reference database. This works best if the reference database is cut to the region of interest. Therefore to run the eDNA expeditions samples, a reference database was built for each target loci, using the [CRABS creating reference databases workflow](https://github.com/gjeunen/reference_database_creator). 

The reference database creator can download sequences directly from NCBI based on a word search, or alternatively can directly download the mitofish database. The workflow consists of 
1. Downloading the ncbi taxonomic database 
2. Downloading the target reference database based on a word search
3. Searching for the target primers in the downloaded sequences
4. Aligning the sequences that had primers with the remaining reference database (because many sequences will not have the primer region available)
5. Building the taxonomic reference file based on the sequences collected

This document shows the searches and analyses made in March 2023 to build the reference databases.

### Set up 

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

### Biomarkers

The different loci analysed in the eDNA expeditions project and the reference databases and searches done for these are:

|Loci |F-primer| R-primer | Reference database | Search text |
|-----|--------|----------|--| ---|
|Mifish-UE| GT**Y**GGTAAA**W**CTCGTGCCAGC      |CATAGTGGGGTATCTAATCC**Y**AGTTTG | NCBI |12S_ncbi_mifish_pga.fasta | 
|Mimammal-UEB  |  GG**RY**TGGT**H**AATTTCGTGCCAGC    | CATAGTG**R**GGTATCTAATC**Y**CAGTTTG | NCBI | 12S_ncbi_mimammals_pga.fasta |
| Teleo - 12S | ACACCGCCCGTCACTCT | CTTCCGGTACACTTACCATG | NCBI | 12S_ncbi_teleo_pga.fasta |
| Leray - COI | GGWACWGGWTGAACWGTWTAYCCYCC | TAIACYTCIGGRTGICCRAARAAYCA | MIDORI | no searches made MIDORI_UNIQ_SP_NUC_GB246_CO1_QIIME.fasta |
| Vert - 16S | AGACGAGAAGACCCYdTGGAGCTT | GATCCAACATCGAGGTCGTAA | NCBI | 16S_ncbi_vert_pga.fasta |

For Mifish and Mimammal the primers used for searching the reference databases and the sequences were a consensus primer from the different versions available. In bold are marked the non-specific basepairs due to building this consensus sequence. 

It is also possible to build the reference database using the mitofish database. In this case we used a search of NCBI.


### Sequence and taxonomy files



The setup steps consist of downloading the databases of choice and the taxonomic dump of ncbi.

#### Taxonomy

```
crabs db_download --source taxonomy

```

#### Sequences

Mifish and teleo (221918 sequences collected)

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

A query for mammals was made separately as a trial though the sequences could have been extracted also from the previous download. A total of 321306 sequences were collected.

```

12S[All Fields] AND Mammalia[Organism] AND ("1"[SLEN] : "50000"[SLEN])  

```
And the same for 16S, the query resulted in about 445298 sequences downloaded.

```

16S[All Fields] AND "Eukaryota"[Organism] AND ("1"[SLEN] : "50000"[SLEN]) AND mitochondrial[All Fields]

```

For COI the following search was made. Originally the MIDORI database was used, but for the sake of comparison, we train a COI database.

```


```

### Extract amplicon sequences

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
| Mifish-UE | ncbi |221918 | 64025 | 85358 |
| Mimammal-UEB | ncbi | 321306 | 140778 | 174006 |
| Teleo | ncbi | (same as mifish) | 62158 | 89796 |
| 16S | ncbi|445298 |  274324 | 322700|
| Mitofish-Mifish | mitofish | 796239  | 16786 | 31441 |

Mitofish added here for comparison, but in the trial data analysis runs, the ncbi extracted data was used. 

### Assign taxonomy

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

### Cleaning of reference sequences

To work with the intended rdp training workflow, the reference sequence collection needs to be cleaned. 
In the original file there are 4443269 sequences.

```
crabs dereplicate --input COI_ncbi_1_50000_pcr_pga_taxon.tsv --output COI_ncbi_1_50000_pcr_pga_taxon_derep.tsv --method uniq_species
```

After dereplication there are 2617162 sequences.

```
crabs seq_cleanup --input COI_ncbi_1_50000_pcr_pga_taxon_derep.tsv --maxns 0 --output COI_ncbi_1_50000_pcr_pga_taxon_derep_clean.tsv --minlen 100 --maxlen 5000 --nans 6 --enviro yes --species yes


filtering parameters:
removing sequences shorter than 100 and longer than 5000
removing sequences containing more than 0 "N"
removing environmental sequences
removing sequences without a species ID
removing sequences with missing information for at least 6 taxonomic levels
found 2617162 number of sequences in ../COI_ncbi_1_50000_pcr_pga_taxon_derep.tsv
removed 1329349 sequences during filtering
written 1287813 sequences to COI_ncbi_1_50000_pcr_pga_taxon_derep_clean.tsv

        found 544 sequences below Minimum count
        found 11 sequences above Maximum count
        found 1329349 sequences with too many Ns
        found 28542 sequences listed as "Environmental"
        found 1195687 sequences listed as sp.
        found 11276 sequences with missing taxonomic information

```

While with keeping environmental sequences, unclassified at species level and those with missing info on 7 taxonomic levels I had: after cleaning there are 2404389 sequences, with 212774 sequences removed. 

Then this final cleaned database can be output in different formats

```
crabs tax_format --input COI_ncbi_1_50000_pcr_pga_taxon_derep_clean.tsv  --output COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp.fasta --format rdp
```


## Building and training the rdp database

The algorithms used are installed with conda rdpTools

```
conda install -c bioconda rdptools
```

Other necessary scripts can be downloaded from here: https://github.com/GLBRC-TeamMicrobiome/python_scripts. 


The first step is to build the taxonomy file so that it has the correct format. This means a Seq-ID (an NCBI id in this case), and seven taxonomic levels, with no empty fields and no convergent names. For example the same genus name cannot exist with two different family names, so a placeholder name must be added to compensate for this. More information in the tutorial of [John Quensen](https://john-quensen.com/tutorials/training-the-rdp-classifier/).

### Clean-up of taxonomic information

So from the file collected previously we separate out the required fields. We take the fasta header, remove the > separate the taxonomy to columns, and remove the "root" field of the taxonomy. 

```
grep "^>" COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp.fasta | sed 's/^>//g' | sed 's/,/\t/'  > COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp.tsv

```

Next we need to replace any empty fields with information. In this case, if there is for example and empty genus, I will add placeholder with g_family, where the 'family' is the name of the family that the genus belongs in. 

```
awk 'BEGIN{OFS=","} {
  split($3, tax, ","),
  for (i = 1, i <= 7, i++) {
    if (tax[i] == "") {
      tax[i] = substr("kpcofgs", i, 1) "_" tax[i - 1],
    }
  }
  print $1, tax[1], tax[2], tax[3], tax[4], tax[5], tax[6], tax[7],
}' COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp.tsv | sed 's/,/\t/g '> COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled.tsv

```
Despite the fact CRABS cleanup was supposed to remove the taxonomies that were fully unknown, there are still many unknown entries in this file. These will be removed

```
grep -v "g_f_o_c_p_k_" COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled.tsv > COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled_nona.tsv

wc -l COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled_nona.tsv
```

This results in 1287579 sequences.  
With less filtering it was 2392878 sequences, from the original 2404388.

Finally, we need to add a header for the taxonomy file to be complete. 

```
sed -i $'1 i\\\nSeq-ID\tSuperkingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled_nona.tsv
```

### Prepare files for training

The files are prepared for training by utilising the downloaded python scripts:

lineage2taxTrain.py and addFullLineage.py

The original scripts are made for python 2. lineage2Train.py gets stuck at lines 34 and 35, which can be commented out. The print function has to be changed to match python 3. After these small changes the script can be run with python 3: 

```
python /home/ubuntu/Saara/github/python_scripts/lineage2taxTrain.py COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled_nona.tsv > ready4train_COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp.tsv

```

Next the modified lineage is added to the sequence file as well.

First remove the previous taxonomy string. This shouldn't matter, but the code fails if this is not removed. 

```
cut -f1 COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp.fasta > COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_names.fasta
```

If there are missing taxonomies, the addFullLineage script will stop. Therefore, next remove all sequences which have already been removed from the taxonomy file, due to no taxonomic string.

```
cut -f1 COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled_nona.tsv > Seq_IDs.tsv

seqtk subseq COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_names.fasta Seq_IDs.tsv  > COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_names2.fasta

```
Make sure that you have the correct number of sequences now

```
grep -c "^>" COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_names2.fasta
1287579

```

Finally, this file can be prepared for training:

```
python /home/ubuntu/Saara/github/python_scripts/addFullLineage.py COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled_nona.tsv COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_names2.fasta > ready4train_COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp.fasta
```

### Train the files

These files are then used for training

```
classifier -Xmx16g train -o training_files -s ready4train_COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp.fasta -t ready4train_COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp.tsv 

```

This gives us an error with a list of convergent names:

```
thalia  genus   2
zonaria genus   2
olea    genus   2
hemineura       genus   2
bostrychia      genus   2
rhizophagus     genus   2
ganonema        genus   2
halopteris      genus   2
turbinaria      genus   2
gustavia        genus   2
drymonia        genus   2
donax   genus   2
trichia genus   2
chilodontidae   family  2
acrotylus       genus   2
solieria        genus   2
lessonia        genus   2
allotropa       genus   2
karenia genus   2
eisenia genus   2
galene  genus   2
adamsiella      genus   2
proteus genus   2
elachista       genus   2
batella genus   2
grania  genus   2
rhodococcus     genus   2
lobophora       genus   2
trichodon       genus   2
metschnikowia   genus   2
ozophora        genus   2
phyllophorella  genus   2
satyrium        genus   2
ptilophora      genus   2
chondria        genus   2
paracoccus      genus   2
byblis  genus   2
heterococcus    genus   2
nemastoma       genus   2
bacillus        genus   2
dracunculus     genus   2
cepheidae       family  2
bartramia       genus   2
darwinella      genus   2
aneura  genus   2
cryptococcus    genus   2
scoparia        genus   2
amesia  genus   2
flammulina      genus   2
agardhiella     genus   2
glaucosphaera   genus   2
cystophora      genus   2
cereus  genus   2
achlya  genus   2
tephrosia       genus   2
xylophagaidae   family  2
chondracanthus  genus   2
heterocheilidae family  2
crassa  genus   2
polypodium      genus   2
ariopsis        genus   2
alaria  genus   2
morus   genus   2
paxillus        genus   2
iris    genus   2
sphaerococcus   genus   2
mastophora      genus   2
dilophus        genus   2
lactarius       genus   2
huonia  genus   2
dacrydium       genus   2
trigonostomum   genus   2

```

To make these genus names unique, I will replace them with a string as such:

f_family_g_genus

For this I use a short script: subset_genus_names.py

First copy the names in a text file, keep only the first column and make the first letter of each genus uppercase.

```
grep "family" genus_names.txt | cut -d " " -f1 | sed 's/.*/\u&/' > family_names_uc.txt
grep "genus" genus_names.txt | cut -d " " -f1 | sed 's/.*/\u&/' > genus_names_uc.txt

```

### Manual check of one outlier (skipped with more filtering)

In this case (before more strong filtering) there is one example where the genus name didn't originally exist, and was replaced by
G_xylophagaidae. This should be removed and replaced manually in the final file.



```
python /home/ubuntu/Saara/github/python_scripts/subset_genus_names.py family_names_uc.txt genus_names_uc.txt COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled_nona.tsv COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled_mod.tsv

```


```
grep "g_Xylophagaidae" COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled_mod.tsv
```

There are 2 orders with the family Xylophagaidae and no further taxonomic information. All have species designations. I will add o_Order_g_family for the replacement string. 

```
OM910837        Eukaryota       Mollusca        Bivalvia        Myida   o_Myida_f_Xylophagaidae g_Xylophagaidae Xylophagaidae_sp._E23
KR390440        Eukaryota       Arthropoda      Insecta Diptera o_Diptera_f_Xylophagaidae       g_Xylophagaidae Xylophagidae_sp._BOLD:ABV1101
```

```
awk 'BEGIN{FS=OFS="\t"} $5=="Myida" && $7=="g_Xylophagaidae" {$7="o_Myida_g_Xylophagaidae"} 1' COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled_mod.tsv > COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled_mod2.tsv

awk 'BEGIN{FS=OFS="\t"} $5=="Diptera" && $7=="g_Xylophagaidae" {$7="o_Diptera_g_Xylophagaidae"} 1' COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled_mod2.tsv > COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled_mod.tsv

```

### Re-run training:

Now prepare the training files and run the training again:

```
python /home/ubuntu/Saara/github/python_scripts/lineage2taxTrain.py COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled_mod.tsv > ready4train_COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp.tsv

```

Just in case, make sure that you have the same sequence names again. Otherwise the script stops. 

```
cut -f1 COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled_mod.tsv > Seq_IDs.tsv

seqtk subseq COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_names.fasta Seq_IDs.tsv  > COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_names2.fasta

```

Then modify the fasta file

```
python /home/ubuntu/Saara/github/python_scripts/addFullLineage.py COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_filled_mod.tsv COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp_names2.fasta > ready4train_COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp.fasta
```


Then train

```
classifier -Xmx40g train -o training_files -s ready4train_COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp.fasta -t ready4train_COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_rdp.tsv 
```

Finally an rRNAClassifier.properties file needs to be copied into the generated folder (training_files). 
This is taken from the terrimporter download for now, and contains the following information:

```
Sample ResourceBundle properties file
bergeyTree=bergeyTrainingTree.xml

probabilityList=genus_wordConditionalProbList.txt

probabilityIndex=wordConditionalProbIndexArr.txt

wordPrior=logWordPrior.txt

classifierVersion=RDP Naive Bayesian rRNA Classifier Version 2.13, August 2021
```

## Testing the trained file with mock data

To test if the trained reference database works, we analyse a mock fish dataset from Hleap et al. 2021. This dataset has been analysed also with vsearch and the terrimporter trained COI dataset. 

The representative sequences can be found in: /home/ubuntu/Saara/github/PacMAN-pipeline/results/Hleap2021/runs/COI_new_pacman/03-dada2/rep-seqs.fna 

```
cd /home/ubuntu/Saara/github/PacMAN-pipeline/results/Hleap2021/runs/COI_new_pacman/03-dada2/

rdp_classifier -Xmx16g classify -t /home/ubuntu/data/databases/COI_ncbi/trial_files/training_files/rRNAClassifier.properties -o COI_ncbi_rdp.output -q rep-seqs.fna

```

The trained dataset was not complete. With the complete dataset the training comes with a memory error until XmX 38g, and stops (crashes?) with 40g. 

For the final analysis we will use the terrimporter trained database. However for the 16S we will need to train our own database, as I have not seen this elsewhere. I will try this, as it is much smaller than the COI database.


## Training a 16S database

To ensure that I have the most up-to-date information I will do the full workflow of downloading the data, pcr, pga, and then formatting for training, (and vsearch).

### Download database

```
cd /home/ubuntu/data/databases/16S/202311
conda activate CRABS

crabs db_download --source ncbi --database nucleotide --query '16S[All Fields] AND "Eukaryota"[Organism] AND ("1"[SLEN] : "50000"[SLEN]) AND mitochondrial[All Fields]' --output 16S_ncbi_euk_1_50000.fasta --keep_original yes --email saara.suominen.work@gmail.com --batchsize 5000

downloading sequences from NCBI
looking up the number of sequences that match the query
found 453175 number of sequences matching the query

found 82 sequences with incorrect accession format
written 453093 sequences to 16S_ncbi_euk_1_50000.fasta

```

### Search and select for primer region

```
crabs insilico_pcr --input 16S_ncbi_euk_1_50000.fasta --output 16S_ncbi_euk_1_50000_pcr.fasta --fwd AGACGAGAAGACCCYdTGGAGCTT --rev GATCCAACATCGAGGTCGTAA --error 4.5

found primers in 264262 sequences
reading sequences without primer-binding regions into memory

reverse complementing 188831 untrimmed sequences
running in silico PCR on 188831 reverse complemented untrimmed sequences
counting the number of sequences found by in silico PCR

found primers in 15260 sequences
```

### Find sequences of the region that have primers cut out

```
crabs pga --input 16S_ncbi_euk_1_50000.fasta --output 16S_ncbi_euk_1_50000_pga.fasta --database 16S_ncbi_euk_1_50000_pcr.fasta --fwd AGACGAGAAGACCCYdTGGAGCTT --rev GATCCAACATCGAGGTCGTAA --speed medium  --percid 0.8 --coverage 0.8 --filter_method relaxed

found 279522 number of sequences in 16S_ncbi_euk_1_50000_pcr.fasta
found 453093 number of sequences in 16S_ncbi_euk_1_50000.fasta

found 169878 (97.87%) number of sequences in 16S_ncbi_euk_1_50000.fasta shorter than 10,000 bp
running pairwise global alignment on 169878 number of sequences. This may take a while!

48841 out of 169878 number of sequences aligned to 16S_ncbi_euk_1_50000_pcr.fasta
48841 number of sequences achieved an alignment that passed filtering thresholds
written 328363 sequences to 16S_ncbi_euk_1_50000_pga.fasta
```

### Add taxonomic information

```
crabs assign_tax --input 16S_ncbi_euk_1_50000_pga.fasta --output 16S_ncbi_euk_1_50000_pga_taxa.tsv --acc2tax /home/ubuntu/data/databases/ncbi/taxonomy20230224/nucl_gb.accession2taxid --taxid /home/ubuntu/data/databases/ncbi/taxonomy20230224/nodes.dmp --name /home/ubuntu/data/databases/ncbi/taxonomy20230224/names.dmp  -w "yes"

found 328363 accession numbers in 16S_ncbi_euk_1_50000_pga.fasta
could not find a taxonomic ID for 1 entries
328362 accession numbers resulted in 88878 unique tax ID numbers
generating taxonomic lineages for 88878 tax ID numbers
assigning a taxonomic lineage to 328363 accession numbers
written 328362 entries to 16S_ncbi_euk_1_50000_pga_taxa.tsv

```

### Clean up sequences

```
crabs dereplicate --input 16S_ncbi_euk_1_50000_pga_taxa.tsv --output 16S_ncbi_euk_1_50000_pga_taxa_derep.tsv --method uniq_species
```
There are 3 possibilities for dereplication 

- strict: only unique sequences will be retained, irrespective of taxonomy
- single_species: for each species in the database, a single sequence is retained
- uniq_species: for each species in the database, all unique sequences are retained


After dereplication there are 155053 sequences.

```
crabs seq_cleanup --input 16S_ncbi_euk_1_50000_pga_taxa_derep.tsv --maxns 0 --output 16S_ncbi_euk_1_50000_pga_taxa_derep_clean.tsv --minlen 100 --maxlen 5000 --nans 6 --enviro yes --species yes

filtering parameters:
removing sequences shorter than 100 and longer than 5000
removing sequences containing more than 0 "N"
removing environmental sequences
removing sequences without a species ID
removing sequences with missing information for at least 6 taxonomic levels
found 155053 number of sequences in 16S_ncbi_euk_1_50000_pga_taxa_derep.tsv
removed 29162 sequences during filtering
written 125891 sequences to 16S_ncbi_euk_1_50000_pga_taxa_derep_clean.tsv

        found 1175 sequences below Minimum count
        found 29162 sequences with too many Ns
        found 163 sequences listed as "Environmental"
        found 23789 sequences listed as sp.
        found 663 sequences with missing taxonomic information

```
 

Then this final cleaned database can be output in different formats

```
crabs tax_format --input 16S_ncbi_euk_1_50000_pga_taxa_derep_clean.tsv  --output 16S_ncbi_euk_1_50000_pga_taxa_derep_clean_rdp.fasta --format rdp

crabs tax_format --input 16S_ncbi_euk_1_50000_pga_taxa_derep_clean.tsv  --output 16S_ncbi_euk_1_50000_pga_taxa_derep_clean_sintax.fasta --format sintax
```

## New 12S databases

There are some 12S databases available also from terrimporter that already trained, however these are very small compared to the amount of information we can collect from ncbi.  

"Created from the MitoFish database [accessed March 16, 2020]. This version contains 2853 reference sequences and 4751 taxa at all ranks."

Or for 12Svert: "This version contains 19,654 reference sequences and 15,007 taxa at all ranks, including 9,564 species"

Therefore we will build also the databases used for the different 12S markers. It is not clear if the database needs to be limited to the region of interest, so we will build it now as such (i.e. including in-silico PCR, and PGA for all markers)

### Download database

Mifish and teleo:

```
./crabs db_download \
--source ncbi \
--database nucleotide \
--query '12S[All Fields] AND mitochondrial[All Fields] AND ("1"[SLEN] : "50000"[SLEN])' \
--output 12S_mito_ncbi_1_50000.fasta \
--keep_original yes \
--email saara.suominen.work@gmail.com \
--batchsize 5000 

downloading sequences from NCBI
looking up the number of sequences that match the query
found 239635 number of sequences matching the query
starting the download

```

Mimammal:

```
crabs db_download \
--source ncbi \
--database nucleotide \
--query '12S[All Fields] AND Mammalia[Organism] AND ("1"[SLEN] : "50000"[SLEN])' \
--output 12S_mammal_ncbi_1_50000.fasta \
--keep_original yes \
--email saara.suominen.work@gmail.com \
--batchsize 5000 

downloading sequences from NCBI
looking up the number of sequences that match the query
found 104287 number of sequences matching the query
```

### Search and select for primer region

Mifish

```
crabs insilico_pcr --input 12S_mito_ncbi_1_50000.fasta --output 12S_mito_ncbi_1_50000_mifish_pcr.fasta --fwd GTYGGTAAAWCTCGTGCCAGC --rev CATAGTGGGGTATCTAATCCYAGTTTG --error 4.5

found primers in 69658 sequences
found primers in 1801 sequences

crabs pga --input 12S_mito_ncbi_1_50000.fasta --output 12S_mito_ncbi_1_50000_mifish_pcr_pga.fasta --database 12S_mito_ncbi_1_50000_mifish_pcr.fasta --fwd GTYGGTAAAWCTCGTGCCAGC --rev CATAGTGGGGTATCTAATCCYAGTTTG --speed medium  --percid 0.8 --coverage 0.8 --filter_method relaxed

25827 number of sequences achieved an alignment that passed filtering thresholds
written 97286 sequences to 12S_mito_ncbi_1_50000_mifish_pcr_pga.fasta
```

Teleo

```
crabs insilico_pcr --input 12S_mito_ncbi_1_50000.fasta --output 12S_mito_ncbi_1_50000_teleo_pcr.fasta --fwd ACACCGCCCGTCACTCT --rev CTTCCGGTACACTTACCATG --error 4.5

found primers in 66642 sequences
found primers in 2593 sequences

crabs pga --input 12S_mito_ncbi_1_50000.fasta --output 12S_mito_ncbi_1_50000_teleo_pcr_pga.fasta --database 12S_mito_ncbi_1_50000_teleo_pcr.fasta --fwd ACACCGCCCGTCACTCT --rev CTTCCGGTACACTTACCATG --speed medium  --percid 0.8 --coverage 0.8 --filter_method relaxed

14699 number of sequences achieved an alignment that passed filtering thresholds
written 83934 sequences to 12S_mito_ncbi_1_50000_teleo_pcr_pga.fasta

```

Mimammal (try with both downloaded databases, most likely there is an overlap)

```
crabs insilico_pcr --input 12S_mito_ncbi_1_50000.fasta --output 12S_mammal_ncbi_1_50000_pcr.fasta --fwd GGRYTGGTHAATTTCGTGCCAGC --rev CATAGTGRGGTATCTAATCYCAGTTTG --error 4.5

found primers in 63064 sequences
found primers in 195 sequences

crabs insilico_pcr --input 12S_mammal_ncbi_1_50000.fasta --output 12S_mammal_ncbi_1_50000_pcr_mammal.fasta --fwd GGRYTGGTHAATTTCGTGCCAGC --rev CATAGTGRGGTATCTAATCYCAGTTTG --error 4.5

found primers in 150985 sequences
found primers in 267 sequences

crabs pga --input 12S_mito_ncbi_1_50000.fasta --output 12S_mammal_ncbi_1_50000_pcr_pga.fasta --database 12S_mammal_ncbi_1_50000_pcr_mammal.fasta --fwd GGRYTGGTHAATTTCGTGCCAGC --rev CATAGTGRGGTATCTAATCYCAGTTTG --speed medium  --percid 0.8 --coverage 0.8 --filter_method relaxed

```

### Add taxonomic information

Mifish

```
crabs assign_tax --input 12S_mito_ncbi_1_50000_mifish_pcr_pga.fasta --output 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa.tsv --acc2tax /home/ubuntu/data/databases/ncbi/taxonomy20230224/nucl_gb.accession2taxid --taxid /home/ubuntu/data/databases/ncbi/taxonomy20230224/nodes.dmp --name /home/ubuntu/data/databases/ncbi/taxonomy20230224/names.dmp  -w "yes"

```

Teleo

```
crabs assign_tax --input 12S_mito_ncbi_1_50000_teleo_pcr_pga.fasta --output 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa.tsv --acc2tax /home/ubuntu/data/databases/ncbi/taxonomy20230224/nucl_gb.accession2taxid --taxid /home/ubuntu/data/databases/ncbi/taxonomy20230224/nodes.dmp --name /home/ubuntu/data/databases/ncbi/taxonomy20230224/names.dmp  -w "yes"

```
Mimammal

```
crabs assign_tax --input 12S_mammal_ncbi_1_50000_pcr_pga.fasta --output 12S_mammal_ncbi_1_50000_pcr_pga_taxa.tsv --acc2tax /home/ubuntu/data/databases/ncbi/taxonomy20230224/nucl_gb.accession2taxid --taxid /home/ubuntu/data/databases/ncbi/taxonomy20230224/nodes.dmp --name /home/ubuntu/data/databases/ncbi/taxonomy20230224/names.dmp  -w "yes"

```



### Clean up sequences

Mifish

```
crabs dereplicate --input 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa.tsv --output 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep.tsv --method uniq_species

found 97286 sequences in 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa.tsv
written 34064 sequences to 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep.tsv

crabs seq_cleanup --input 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep.tsv --maxns 0 --output 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean.tsv --minlen 100 --maxlen 5000 --nans 6 --enviro yes --species yes

written 29667 sequences to 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean.tsv

crabs tax_format --input 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean.tsv  --output 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp.fasta --format rdp

crabs tax_format --input 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean.tsv --output 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_sintax.fasta --format sintax
```

Teleo

```
crabs dereplicate --input 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa.tsv --output 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_derep.tsv --method uniq_species

found 83934 sequences in 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa.tsv
written 27906 sequences to 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_derep.tsv

crabs seq_cleanup --input 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_derep.tsv --maxns 0 --output 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_derep_clean.tsv --minlen 50 --maxlen 5000 --nans 6 --enviro yes --species yes

written 21789 sequences to 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_derep_clean.tsv

crabs tax_format --input 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_derep_clean.tsv  --output 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_clean_rdp.fasta --format rdp

crabs tax_format --input 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_derep_clean.tsv  --output 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_clean_sintax.fasta --format sintax

```


Mimammal

```
crabs dereplicate --input 12S_mammal_ncbi_1_50000_pcr_pga_taxa.tsv --output 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep.tsv --method uniq_species

found 148195 sequences in 12S_mammal_ncbi_1_50000_pcr_pga_taxa.tsv
written 35445 sequences to 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep.tsv

crabs seq_cleanup --input 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep.tsv --maxns 0 --output 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean.tsv --minlen 100 --maxlen 5000 --nans 6 --enviro yes --species yes

written 30520 sequences to 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean.tsv

crabs tax_format --input 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean.tsv  --output 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp.fasta --format rdp

crabs tax_format --input 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean.tsv --output 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_sintax.fasta --format sintax

```


### Clean-up of taxonomic information

Mifish

```
grep "^>" 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp.fasta | sed 's/^>//g' | sed 's/,/\t/'  > 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp.tsv

awk 'BEGIN{OFS=","} {
  split($3, tax, ","),
  for (i = 1, i <= 7, i++) {
    if (tax[i] == "") {
      tax[i] = substr("kpcofgs", i, 1) "_" tax[i - 1],
    }
  }
  print $1, tax[1], tax[2], tax[3], tax[4], tax[5], tax[6], tax[7],
}' 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp.tsv | sed 's/,/\t/g '> 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp_filled.tsv

grep -v "g_f_o_c_p_k_" 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp_filled.tsv > 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp_filled_nona.tsv

wc -l 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp_filled_nona.tsv

sed -i $'1 i\\\nSeq-ID\tSuperkingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp_filled_nona.tsv
```

Teleo

```
grep "^>" 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_clean_rdp.fasta | sed 's/^>//g' | sed 's/,/\t/'  > 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_clean_rdp.tsv

awk 'BEGIN{OFS=","} {
  split($3, tax, ","),
  for (i = 1, i <= 7, i++) {
    if (tax[i] == "") {
      tax[i] = substr("kpcofgs", i, 1) "_" tax[i - 1],
    }
  }
  print $1, tax[1], tax[2], tax[3], tax[4], tax[5], tax[6], tax[7],
}' 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_clean_rdp.tsv | sed 's/,/\t/g '> 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_clean_rdp_filled.tsv

grep -v "g_f_o_c_p_k_" 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_clean_rdp_filled.tsv > 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_clean_rdp_filled_nona.tsv

wc -l 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_clean_rdp_filled_nona.tsv

sed -i $'1 i\\\nSeq-ID\tSuperkingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_clean_rdp_filled_nona.tsv
```

Mimammal

```
grep "^>" 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp.fasta | sed 's/^>//g' | sed 's/,/\t/'  > 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp.tsv

awk 'BEGIN{OFS=","} {
  split($3, tax, ","),
  for (i = 1, i <= 7, i++) {
    if (tax[i] == "") {
      tax[i] = substr("kpcofgs", i, 1) "_" tax[i - 1],
    }
  }
  print $1, tax[1], tax[2], tax[3], tax[4], tax[5], tax[6], tax[7],
}' 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp.tsv | sed 's/,/\t/g '> 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp_filled.tsv

grep -v "g_f_o_c_p_k_" 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp_filled.tsv > 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp_filled_nona.tsv

wc -l 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp_filled_nona.tsv

sed -i $'1 i\\\nSeq-ID\tSuperkingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp_filled_nona.tsv
```


### Prepare files for training

Mifish

```
python /home/ubuntu/Saara/github/python_scripts/lineage2taxTrain.py 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp_filled_nona.tsv > ready4train_12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp.tsv

cut -f1 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp.fasta > 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp_names.fasta

cut -f1 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp_filled_nona.tsv > Seq_IDs.tsv

seqtk subseq 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp_names.fasta Seq_IDs.tsv  > 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp_names2.fasta

grep -c "^>" 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp_names2.fasta
29665

python /home/ubuntu/Saara/github/python_scripts/addFullLineage.py 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp_filled_nona.tsv 12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp_names2.fasta > ready4train_12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp.fasta
```

Train the files:

```
classifier -Xmx16g train -o training_files_mifish -s ready4train_12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp.fasta -t ready4train_12S_mito_ncbi_1_50000_mifish_pcr_pga_taxa_derep_clean_rdp.tsv

```


Teleo

```
python /home/ubuntu/Saara/github/python_scripts/lineage2taxTrain.py 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_clean_rdp_filled_nona.tsv > ready4train_12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_derep_clean_rdp.tsv

cut -f1 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_clean_rdp.fasta  > 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_derep_clean_rdp_names.fasta

cut -f1 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_clean_rdp_filled_nona.tsv > Seq_IDs.tsv

seqtk subseq 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_derep_clean_rdp_names.fasta Seq_IDs.tsv  > 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_derep_clean_rdp_names2.fasta

grep -c "^>" 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_derep_clean_rdp_names2.fasta
21788

python /home/ubuntu/Saara/github/python_scripts/addFullLineage.py 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_clean_rdp_filled_nona.tsv 12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_derep_clean_rdp_names2.fasta > ready4train_12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_derep_clean_rdp.fasta
```

Train the files:

```
classifier -Xmx16g train -o training_files_teleo -s ready4train_12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_derep_clean_rdp.fasta -t ready4train_12S_mito_ncbi_1_50000_teleo_pcr_pga_taxa_derep_clean_rdp.tsv

```

Mimammal

```
python /home/ubuntu/Saara/github/python_scripts/lineage2taxTrain.py 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp_filled_nona.tsv > ready4train_12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp.tsv

cut -f1 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp.fasta > 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp_names.fasta

cut -f1 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp_filled_nona.tsv > Seq_IDs.tsv

seqtk subseq 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp_names.fasta Seq_IDs.tsv  > 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp_names2.fasta

grep -c "^>" 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp_names2.fasta
30519

python /home/ubuntu/Saara/github/python_scripts/addFullLineage.py 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp_filled_nona.tsv 12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp_names2.fasta > ready4train_12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp.fasta
```

Train the files:

```
classifier -Xmx16g train -o training_files_mimammal -s ready4train_12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp.fasta -t ready4train_12S_mammal_ncbi_1_50000_pcr_pga_taxa_derep_taxa_clean_rdp.tsv

```