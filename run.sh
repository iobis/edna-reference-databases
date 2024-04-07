#!/bin/bash

primer_12s_mifish_f="GTYGGTAAAWCTCGTGCCAGC"
primer_12s_mifish_r="CATAGTGGGGTATCTAATCCYAGTTTG"
primer_12s_mimammal_f="GGRYTGGTHAATTTCGTGCCAGC"
primer_12s_mimammal_r="CATAGTGRGGTATCTAATCYCAGTTTG"
primer_12s_teleo_f="ACACCGCCCGTCACTCT"
primer_12s_teleo_r="CTTCCGGTACACTTACCATG"
primer_coi_f="GGWACWGGWTGAACWGTWTAYCCYCC"
primer_coi_r="TAIACYTCIGGRTGICCRAARAAYCA"
primer_16s_f="AGACGAGAAGACCCYdTGGAGCTT"
primer_16s_r="GATCCAACATCGAGGTCGTAA"

threads=6

if ! [ -f nucl_gb.accession2taxid ]; then
crabs db_download --source taxonomy
fi

if ! [ -f coi_ncbi_1_50000.fasta ]; then
crabs db_download --source ncbi --database nucleotide --query 'COI[All Fields] OR CO1[All Fields] OR cytochrome oxidase subunit I[All Fields] OR cytochrome oxidase subunit 1[All Fields] OR cytochrome c oxidase subunit I[All Fields] OR cytochrome c oxidase subunit 1[All Fields] OR COX1[All Fields] AND ("50"[SLEN] : "50000"[SLEN])' --output coi_ncbi_1_50000.fasta --keep_original yes --email s@gmail.com --batchsize 5000
crabs db_download --source ncbi --database nucleotide --query '12S[All Fields] AND mitochondrial[All Fields] AND ("1"[SLEN] : "50000"[SLEN])' --output 12s_ncbi_1_50000.fasta --keep_original yes --email s@gmail.com --batchsize 5000
crabs db_download --source ncbi --database nucleotide --query '16S[All Fields] AND "Eukaryota"[Organism] AND ("1"[SLEN] : "50000"[SLEN]) AND mitochondrial[All Fields]' --output 16s_ncbi_1_50000.fasta --keep_original yes --email s@gmail.com --batchsize 5000
fi

if ! [ -f coi_ncbi.fasta ]; then
crabs insilico_pcr --threads "$threads" --input coi_ncbi_1_50000.fasta --output coi_ncbi.fasta --fwd "$primer_coi_f" --rev "$primer_coi_r" --error 4.5
crabs insilico_pcr --threads "$threads" --input 12s_ncbi_1_50000.fasta --output 12s_ncbi_mifish.fasta --fwd "$primer_12s_mifish_f" --rev "$primer_12s_mifish_r" --error 4.5
crabs insilico_pcr --threads "$threads" --input 12s_ncbi_1_50000.fasta --output 12s_ncbi_mimammal.fasta --fwd "$primer_12s_mimammal_f" --rev "$primer_12s_mimammal_r" --error 4.5
crabs insilico_pcr --threads "$threads" --input 12s_ncbi_1_50000.fasta --output 12s_ncbi_teleo.fasta --fwd "$primer_12s_teleo_f" --rev "$primer_12s_teleo_r" --error 4.5
crabs insilico_pcr --threads "$threads" --input 16s_ncbi_1_50000.fasta --output 16s_ncbi.fasta --fwd "$primer_16s_f" --rev "$primer_16s_r" --error 4.5
fi

if ! [ -f coi_ncbi_pga.fasta ]; then
crabs pga --input coi_ncbi_1_50000.fasta --output coi_ncbi_pga.fasta --database coi_ncbi.fasta --fwd "$primer_coi_f" --rev "$primer_coi_r" --speed medium --percid 0.8 --coverage 0.8 --filter_method relaxed
crabs pga --input 12s_ncbi_1_50000.fasta --output 12s_ncbi_mifish_pga.fasta --database 12s_ncbi_mifish.fasta --fwd "$primer_12s_mifish_f" --rev "$primer_12s_mifish_r" --speed medium --percid 0.8 --coverage 0.8 --filter_method relaxed
crabs pga --input 12s_ncbi_1_50000.fasta --output 12s_ncbi_mimammal_pga.fasta --database 12s_ncbi_mimammal.fasta --fwd "$primer_12s_mimammal_f" --rev "$primer_12s_mimammal_r" --speed medium --percid 0.8 --coverage 0.8 --filter_method relaxed
crabs pga --input 12s_ncbi_1_50000.fasta --output 12s_ncbi_teleo_pga.fasta --database 12s_ncbi_teleo.fasta --fwd "$primer_12s_teleo_f" --rev "$primer_12s_teleo_r" --speed medium --percid 0.8 --coverage 0.8 --filter_method relaxed
crabs pga --input 16s_ncbi_1_50000.fasta --output 16s_ncbi_pga.fasta --database 16s_ncbi.fasta --fwd "$primer_16s_f" --rev "$primer_16s_r" --speed medium --percid 0.8 --coverage 0.8 --filter_method relaxed
fi

if ! [ -f coi_ncbi_pga_taxa.tsv ]; then
crabs assign_tax --input coi_ncbi_pga.fasta --output coi_ncbi_pga_taxa.tsv --acc2tax nucl_gb.accession2taxid  --taxid nodes.dmp --name names.dmp
crabs assign_tax --input 12s_ncbi_mifish_pga.fasta --output 12s_ncbi_mifish_pga_taxa.tsv --acc2tax nucl_gb.accession2taxid  --taxid nodes.dmp --name names.dmp
crabs assign_tax --input 12s_ncbi_mimammal_pga.fasta --output 12s_ncbi_mimammal_pga_taxa.tsv --acc2tax nucl_gb.accession2taxid  --taxid nodes.dmp --name names.dmp
crabs assign_tax --input 12s_ncbi_teleo_pga.fasta --output 12s_ncbi_teleo_pga_taxa.tsv --acc2tax nucl_gb.accession2taxid  --taxid nodes.dmp --name names.dmp
crabs assign_tax --input 16s_ncbi_pga.fasta --output 16s_ncbi_pga_taxa.tsv --acc2tax nucl_gb.accession2taxid  --taxid nodes.dmp --name names.dmp
fi

if ! [ -f coi_ncbi_pga_pacman.tsv ]; then
awk -F'\t' '{OFS="\t"; $10=$3";"$4";"$5";"$6";"$7";"$8";"$9; print}' coi_ncbi_pga_taxa.tsv > coi_ncbi_pga_pacman.tsv
awk -F'\t' '{OFS="\t"; $10=$3";"$4";"$5";"$6";"$7";"$8";"$9; print}' 12s_ncbi_mifish_pga_taxa.tsv > 12s_ncbi_mifish_pga_pacman.tsv
awk -F'\t' '{OFS="\t"; $10=$3";"$4";"$5";"$6";"$7";"$8";"$9; print}' 12s_ncbi_mimammal_pga_taxa.tsv > 12s_ncbi_mimammal_pga_pacman.tsv
awk -F'\t' '{OFS="\t"; $10=$3";"$4";"$5";"$6";"$7";"$8";"$9; print}' 12s_ncbi_teleo_pga_taxa.tsv > 12s_ncbi_teleo_pga_pacman.tsv
awk -F'\t' '{OFS="\t"; $10=$3";"$4";"$5";"$6";"$7";"$8";"$9; print}' 16s_ncbi_pga_taxa.tsv > 16s_ncbi_pga_pacman.tsv
fi

if ! [ -f coi_ncbi_pga_taxa_pacmanformat.tsv ]; then
cut -f1,10 coi_ncbi_pga_pacman.tsv > coi_ncbi_pga_taxa_pacmanformat.tsv
cut -f1,10 12s_ncbi_mifish_pga_pacman.tsv > 12s_ncbi_mifish_pga_taxa_pacmanformat.tsv
cut -f1,10 12s_ncbi_mimammal_pga_pacman.tsv > 12s_ncbi_mimammal_pga_taxa_pacmanformat.tsv
cut -f1,10 12s_ncbi_teleo_pga_pacman.tsv > 12s_ncbi_teleo_pga_taxa_pacmanformat.tsv
cut -f1,10 16s_ncbi_pga_pacman.tsv > 16s_ncbi_pga_taxa_pacmanformat.tsv
fi

if ! [ -f coi_ncbi_pga_taxa_derep.tsv ]; then
crabs dereplicate --input coi_ncbi_pga_taxa.tsv --output coi_ncbi_pga_taxa_derep.tsv --method uniq_species
crabs dereplicate --input 12s_ncbi_mifish_pga_taxa.tsv --output 12s_ncbi_mifish_pga_taxa_derep.tsv --method uniq_species
crabs dereplicate --input 12s_ncbi_mimammal_pga_taxa.tsv --output 12s_ncbi_mimammal_pga_taxa_derep.tsv --method uniq_species
crabs dereplicate --input 12s_ncbi_teleo_pga_taxa.tsv --output 12s_ncbi_teleo_pga_taxa_derep.tsv --method uniq_species
crabs dereplicate --input 16s_ncbi_pga_taxa.tsv --output 16s_ncbi_pga_taxa_derep.tsv --method uniq_species
fi

if ! [ -f coi_ncbi_pga_taxa_derep_clean.tsv ]; then
crabs seq_cleanup --input coi_ncbi_pga_taxa_derep.tsv --maxns 0 --output coi_ncbi_pga_taxa_derep_clean.tsv --minlen 100 --maxlen 5000 --nans 6 --enviro yes --species yes
crabs seq_cleanup --input 12s_ncbi_mifish_pga_taxa_derep.tsv --maxns 0 --output 12s_ncbi_mifish_pga_taxa_derep_clean.tsv --minlen 100 --maxlen 5000 --nans 6 --enviro yes --species yes
crabs seq_cleanup --input 12s_ncbi_mimammal_pga_taxa_derep.tsv --maxns 0 --output 12s_ncbi_mimammal_pga_taxa_derep_clean.tsv --minlen 100 --maxlen 5000 --nans 6 --enviro yes --species yes
crabs seq_cleanup --input 12s_ncbi_teleo_pga_taxa_derep.tsv --maxns 0 --output 12s_ncbi_teleo_pga_taxa_derep_clean.tsv --minlen 100 --maxlen 5000 --nans 6 --enviro yes --species yes
crabs seq_cleanup --input 16s_ncbi_pga_taxa_derep.tsv --maxns 0 --output 16s_ncbi_pga_taxa_derep_clean.tsv --minlen 100 --maxlen 5000 --nans 6 --enviro yes --species yes
fi

if ! [ -f coi_ncbi_pga_taxa_derep_clean_rdp.fasta ]; then
crabs tax_format --input coi_ncbi_pga_taxa_derep_clean.tsv  --output coi_ncbi_pga_taxa_derep_clean_rdp.fasta --format rdp
crabs tax_format --input 12s_ncbi_mifish_pga_taxa_derep_clean.tsv  --output 12s_ncbi_mifish_pga_taxa_derep_clean_rdp.fasta --format rdp
crabs tax_format --input 12s_ncbi_mimammal_pga_taxa_derep_clean.tsv  --output 12s_ncbi_mimammal_pga_taxa_derep_clean_rdp.fasta --format rdp
crabs tax_format --input 12s_ncbi_teleo_pga_taxa_derep_clean.tsv  --output 12s_ncbi_teleo_pga_taxa_derep_clean_rdp.fasta --format rdp
crabs tax_format --input 16s_ncbi_pga_taxa_derep_clean.tsv  --output 16s_ncbi_pga_taxa_derep_clean_rdp.fasta --format rdp
fi
