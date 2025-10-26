#!/bin/bash

# Input files
annotation_file="/BioII/lulab_b/huangkeyun/zhangys/RNA_locator/references/exo_circRNAs_anno.txt"
genome_file="/BioII/lulab_b/huangkeyun/zhangys/RNA_locator/references/GRCh38.p14.genome.fa"
output_file="./exo_circRNA_sequences.csv"

# Extract BED file (including strand information)
awk -F'\t' 'NR>1 {split($3, pos, "[:-]"); print pos[1] "\t" pos[2]-1 "\t" pos[3] "\t" $1 "\t" "." "\t" $4}' "$annotation_file" > circRNA.bed

# Use bedtools to extract sequences (considering positive and negative strands)
bedtools getfasta -fi "$genome_file" -bed circRNA.bed -fo circRNA_sequences.fa -name -s

# Convert Fasta to CSV
awk 'BEGIN {OFS=","} NR%2==1 {id=substr($0, 2)} NR%2==0 {print id, $0}' circRNA_sequences.fa > "$output_file"

echo "circRNA sequence extraction completed, results saved to $output_file"
