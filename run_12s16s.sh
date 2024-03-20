#!/bin/bash

# This creates 12S reference databases based on a combined 12S and 16S download from NCBI.

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

threads=4

# crabs db_download --source taxonomy

# crabs db_download --source ncbi --database nucleotide --query '12S[All Fields] AND ribosomal[All Fields]  AND ("1"[SLEN] : "5000"[SLEN])' --output 12s_ncbi_1_5000.fasta --keep_original yes --email s@gmail.com --batchsize 5000
# crabs db_download --source ncbi --database nucleotide --query '16S[All Fields] AND ribosomal[All Fields]  AND ("1"[SLEN] : "5000"[SLEN])' --output 16s_ncbi_1_5000.fasta --keep_original yes --email s@gmail.com --batchsize 5000

# cat 12s_ncbi_1_5000.fasta > 12s_16s_ncbi_1_5000_combined.fasta
# echo "" >> 12s_16s_ncbi_1_5000_combined.fasta
# cat 16s_ncbi_1_5000.fasta >> 12s_16s_ncbi_1_5000_combined.fasta

# seqkit rmdup -n < 12s_16s_ncbi_1_5000_combined.fasta > 12s_16s_ncbi_1_5000.fasta

# crabs insilico_pcr --threads "$threads" --input 12s_16s_ncbi_1_5000.fasta --output 12s_16s_ncbi_mifish.fasta --fwd "$primer_12s_mifish_f" --rev "$primer_12s_mifish_r" --error 4.5
# crabs insilico_pcr --threads "$threads" --input 12s_16s_ncbi_1_5000.fasta --output 12s_16s_ncbi_mimammal.fasta --fwd "$primer_12s_mimammal_f" --rev "$primer_12s_mimammal_r" --error 4.5
# crabs insilico_pcr --threads "$threads" --input 12s_16s_ncbi_1_5000.fasta --output 12s_16s_ncbi_teleo.fasta --fwd "$primer_12s_teleo_f" --rev "$primer_12s_teleo_r" --error 4.5

# crabs pga --input 12s_16s_ncbi_1_5000.fasta --output 12s_16s_ncbi_mifish_pga.fasta --database 12s_16s_ncbi_mifish.fasta --fwd "$primer_12s_mifish_f" --rev "$primer_12s_mifish_r" --speed medium --percid 0.8 --coverage 0.8 --filter_method relaxed
# crabs pga --input 12s_16s_ncbi_1_5000.fasta --output 12s_16s_ncbi_mimammal_pga.fasta --database 12s_16s_ncbi_mimammal.fasta --fwd "$primer_12s_mimammal_f" --rev "$primer_12s_mimammal_r" --speed medium --percid 0.8 --coverage 0.8 --filter_method relaxed
# crabs pga --input 12s_16s_ncbi_1_5000.fasta --output 12s_16s_ncbi_teleo_pga.fasta --database 12s_16s_ncbi_teleo.fasta --fwd "$primer_12s_teleo_f" --rev "$primer_12s_teleo_r" --speed medium --percid 0.8 --coverage 0.8 --filter_method relaxed

# crabs assign_tax --input 12s_16s_ncbi_mifish_pga.fasta --output 12s_16s_ncbi_mifish_pga_taxa.tsv --missing 12s_16s_ncbi_mifish_pga_missing_taxa.tsv --acc2tax nucl_gb.accession2taxid  --taxid nodes.dmp --name names.dmp
# crabs assign_tax --input 12s_16s_ncbi_mimammal_pga.fasta --output 12s_16s_ncbi_mimammal_pga_taxa.tsv --missing 12s_16s_ncbi_mimammal_pga_missing_taxa.tsv --acc2tax nucl_gb.accession2taxid  --taxid nodes.dmp --name names.dmp
# crabs assign_tax --input 12s_16s_ncbi_teleo_pga.fasta --output 12s_16s_ncbi_teleo_pga_taxa.tsv --missing 12s_16s_ncbi_teleo_pga_missing_taxa.tsv --acc2tax nucl_gb.accession2taxid  --taxid nodes.dmp --name names.dmp
