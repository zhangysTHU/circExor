#!/bin/bash

if [ ! -d "./Motif_Comparison_Output" ]; then
  mkdir Motif_Comparison_Output
fi

seqkit seq  --dna2rna ./consensus_kmers/EV/motifs.fa > ./Motif_Comparison_Output/EV_motif_consensus_seq_rna.fa
rna2meme ./Motif_Comparison_Output/EV_motif_consensus_seq_rna.fa > ./Motif_Comparison_Output/EV_motif_consensus_seq_rna.meme

tomtom -no-ssc -oc . -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 0.1 \
-norc ./Motif_Comparison_Output/EV_motif_consensus_seq_rna.meme ./motif_databases/CISBP-RNA/Homo_sapiens.meme -o ./Motif_Comparison_Output/homo_EV_tomtom

tomtom -no-ssc -oc . -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 0.1 \
-norc ./Motif_Comparison_Output/cyto_motif_consensus_seq_rna.meme ./motif_databases/CISBP-RNA/Mus_musculus.meme -o ./Motif_Comparison_Output/mus_EV_tomtom

# Remove redundant and keep the best match
grep "M*_0.6" ./Motif_Comparison_Output/EV_tomtom/tomtom.tsv  | sort -k2,2 -k6,6n -u | awk '!a[$2]++'  > ./Motif_Comparison_Output/EV_favour_motif_tomtom.tsv


RBP="./motif_databases/CISBP-RNA/Homo_sapiens.meme"

# get EV_Related RBP name
touch EV_Related_RBP.temp
cat ./Motif_Comparison_Output/EV_favour_motif_tomtom.tsv | grep "M*_0.6" | while read line 
do
    Target_ID=`echo "$line" | cut -f2`
    grep -w $Target_ID $RBP | cut -d" " -f3 >>  EV_Related_RBP.temp
done
paste EV_Related_RBP.temp ./Motif_Comparison_Output/EV_favour_motif_tomtom.tsv > ./Motif_Comparison_Output/EV_favour_motif_final.tsv #


rm EV_Related_RBP.temp 