#!/bin/bash
mkdir -p ./assmbled_kmers
mkdir -p ./consensus_kmers/cyto
mkdir -p ./consensus_kmers/EV
mkdir -p ./circRNA_ligated

flank=10
max_kmer_len=5
# Calculate the length to be appended
add_len=$((flank + ( (max_kmer_len + 1) / 2 )))  # Integer division in bash, use (max_kmer_len+1)/2 for rounding up

# Concatenate the head and tail of circRNA sequences
out_dir=./circRNA_ligated
inputs=(
"../sample_preprocessing/circRNA/Cyto_sequences.fasta"
"../sample_preprocessing/circRNA/EV_sequences.fasta"
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

# Assemble kmers based on kmer and circRNA sequences
python ./motif_analysis/kmer_assmble.py --kmers ../SHAP/cyto_kmers.fasta --circ ./motif_analysis/circRNA_ligated/Cyto_sequences.fasta --out ./assmbled_kmers/cyto_nt.fasta --dedup
python ./motif_analysis/kmer_assmble.py --kmers ../SHAP/EV_kmers.fasta --circ ./motif_analysis/circRNA_ligated/EV_sequences.fasta --out ./assmbled_kmers/EV_nt.fasta --dedup

# Perform multiple sequence alignment on the assembled kmers and generate consensus sequences
python ./motif_analysis/MSA_consensus.py -i ./assmbled_kmers/cyto_nt.fasta -o ./consensus_kmers/cyto
python ./motif_analysis/MSA_consensus.py -i ./assmbled_kmers/EV_nt.fasta -o ./consensus_kmers/EV

# Compare the consensus sequences to known RBP motifs using TomTom
bash ./motif_analysis/cyto_motif_comparision.sh
bash ./motif_analysis/EV_motif_comparison.sh