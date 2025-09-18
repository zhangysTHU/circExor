#!/bin/bash
mkdir -p ./assmbled_kmers
mkdir -p ./consensus_kmers/cyto
mkdir -p ./consensus_kmers/EV
mkdir -p ./circRNA_ligated

flank=10
max_kmer_len=5
# 计算需要追加的长度
add_len=$((flank + ( (max_kmer_len + 1) / 2 )))  # bash 中整数除法，向上取整用 (max_kmer_len+1)/2

# 对 circRNA 序列进行头尾拼接
out_dir=./circRNA_ligated
inputs=(
"/lulabdata3/huangkeyun/zhangys/RNA_locator/ML_python_scripts/reference_preprocessing/circRNA/Cyto_sequences.fasta"
"/lulabdata3/huangkeyun/zhangys/RNA_locator/ML_python_scripts/reference_preprocessing/circRNA/EV_sequences.fasta"
)
for infile in "${inputs[@]}"; do
    fname=$(basename $infile)
    outfile="$out_dir/$fname"

    awk -v add_len=$add_len '
    BEGIN{FS=""; OFS=""}
    /^>/{print; next}
    {
        seq=$0
        seqlen=length(seq)
        if(seqlen >= add_len){
            extra=substr(seq, 1, add_len)
        } else {
            repeat_times=int(add_len / seqlen)
            remainder=add_len % seqlen
            extra=""
            for(i=0;i<repeat_times;i++){extra=extra seq}
            if(remainder>0){extra=extra substr(seq,1,remainder)}
        }
        print seq extra
    }
    ' $infile > $outfile
done
echo "circRNA sequences ligated and saved to $out_dir"

# 根据kmer和circRNA序列进行kmer组装
python /lulabdata3/huangkeyun/zhangys/RNA_locator/motif_analysis/kmer_assmble.py --kmers /lulabdata3/huangkeyun/zhangys/RNA_locator/ML_python_scripts/SHAP/cyto_kmers.fasta --circ /lulabdata3/huangkeyun/zhangys/RNA_locator/motif_analysis/circRNA_ligated/Cyto_sequences.fasta --out ./assmbled_kmers/cyto_nt.fasta --dedup
python /lulabdata3/huangkeyun/zhangys/RNA_locator/motif_analysis/kmer_assmble.py --kmers /lulabdata3/huangkeyun/zhangys/RNA_locator/ML_python_scripts/SHAP/EV_kmers.fasta --circ /lulabdata3/huangkeyun/zhangys/RNA_locator/motif_analysis/circRNA_ligated/EV_sequences.fasta --out ./assmbled_kmers/EV_nt.fasta --dedup

# 对组装后的kmer进行多序列比对并生成共识序列
python /lulabdata3/huangkeyun/zhangys/RNA_locator/motif_analysis/MSA_consensus.py -i ./assmbled_kmers/cyto_nt.fasta -o ./consensus_kmers/cyto
python /lulabdata3/huangkeyun/zhangys/RNA_locator/motif_analysis/MSA_consensus.py -i ./assmbled_kmers/EV_nt.fasta -o ./consensus_kmers/EV

# 将公示序列通过TomTom比对到已知RBP motif
bash /lulabdata3/huangkeyun/zhangys/RNA_locator/motif_analysis/cyto_motif_comparision.sh
bash /lulabdata3/huangkeyun/zhangys/RNA_locator/motif_analysis/EV_motif_comparison.sh