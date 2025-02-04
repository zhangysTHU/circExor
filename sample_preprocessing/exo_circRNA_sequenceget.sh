#!/bin/bash

# 输入文件
annotation_file="/BioII/lulab_b/huangkeyun/zhangys/RNA_locator/references/exo_circRNAs_anno.txt"
genome_file="/BioII/lulab_b/huangkeyun/zhangys/RNA_locator/references/GRCh38.p14.genome.fa"
output_file="./exo_circRNA_sequences.csv"

# 提取 BED 文件（含 strand 信息）
awk -F'\t' 'NR>1 {split($3, pos, "[:-]"); print pos[1] "\t" pos[2]-1 "\t" pos[3] "\t" $1 "\t" "." "\t" $4}' "$annotation_file" > circRNA.bed

# 使用 bedtools 提取序列（考虑正负链）
bedtools getfasta -fi "$genome_file" -bed circRNA.bed -fo circRNA_sequences.fa -name -s

# 转换 Fasta 为 CSV
awk 'BEGIN {OFS=","} NR%2==1 {id=substr($0, 2)} NR%2==0 {print id, $0}' circRNA_sequences.fa > "$output_file"

echo "circRNA 序列提取完成，结果保存至 $output_file"
