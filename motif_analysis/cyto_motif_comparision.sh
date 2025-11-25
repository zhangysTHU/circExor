#!/bin/bash

# Create main output directory
if [ ! -d "./motif_analysis/Motif_Comparison_Output" ]; then
  mkdir Motif_Comparison_Output
fi

# Create species-specific output directories
mkdir -p ./motif_analysis/Motif_Comparison_Output/homo_final
mkdir -p ./motif_analysis/Motif_Comparison_Output/mus_final

seqkit seq  --dna2rna ./motif_analysis/consensus_kmers/cyto/motifs.fa > ./motif_analysis/Motif_Comparison_Output/cyto_motif_consensus_seq_rna.fa
rna2meme ./motif_analysis/Motif_Comparison_Output/cyto_motif_consensus_seq_rna.fa > ./motif_analysis/Motif_Comparison_Output/cyto_motif_consensus_seq_rna.meme

tomtom -no-ssc -oc . -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 0.1 \
-norc ./motif_analysis/Motif_Comparison_Output/cyto_motif_consensus_seq_rna.meme ./motif_analysis/motif_databases/CISBP-RNA/Homo_sapiens.meme -o ./motif_analysis/Motif_Comparison_Output/homo_cyto_tomtom

tomtom -no-ssc -oc . -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 0.1 \
-norc ./motif_analysis/Motif_Comparison_Output/cyto_motif_consensus_seq_rna.meme ./motif_analysis/motif_databases/CISBP-RNA/Mus_musculus.meme -o ./motif_analysis/Motif_Comparison_Output/mus_cyto_tomtom

# Process human results
grep "M*_0.6" ./motif_analysis/Motif_Comparison_Output/homo_cyto_tomtom/tomtom.tsv | sort -k2,2 -k6,6n -u | awk '!a[$2]++'  > ./motif_analysis/Motif_Comparison_Output/homo_final/homo_cyto_favour_motif_tomtom.tsv

# Get human RBP names
touch ./motif_analysis/Motif_Comparison_Output/homo_final/homo_cyto_Related_RBP.temp
cat ./motif_analysis/Motif_Comparison_Output/homo_final/homo_cyto_favour_motif_tomtom.tsv | grep "M*_0.6" | while read line 
do
    Target_ID=`echo "$line" | cut -f2`
    grep -w $Target_ID ./motif_analysis/motif_databases/CISBP-RNA/Homo_sapiens.meme | cut -d" " -f3 >> ./motif_analysis/Motif_Comparison_Output/homo_final/homo_cyto_Related_RBP.temp
done
paste ./motif_analysis/Motif_Comparison_Output/homo_final/homo_cyto_Related_RBP.temp ./motif_analysis/Motif_Comparison_Output/homo_final/homo_cyto_favour_motif_tomtom.tsv > ./motif_analysis/Motif_Comparison_Output/homo_final/homo_cyto_favour_motif_final.tsv
rm ./motif_analysis/Motif_Comparison_Output/homo_final/homo_cyto_Related_RBP.temp

# Process mouse results
grep "M*_0.6" ./motif_analysis/Motif_Comparison_Output/mus_cyto_tomtom/tomtom.tsv | sort -k2,2 -k6,6n -u | awk '!a[$2]++'  > ./motif_analysis/Motif_Comparison_Output/mus_final/mus_cyto_favour_motif_tomtom.tsv

# Get mouse RBP names
touch ./motif_analysis/Motif_Comparison_Output/mus_final/mus_cyto_Related_RBP.temp
cat ./motif_analysis/Motif_Comparison_Output/mus_final/mus_cyto_favour_motif_tomtom.tsv | grep "M*_0.6" | while read line 
do
    Target_ID=`echo "$line" | cut -f2`
    grep -w $Target_ID ./motif_analysis/motif_databases/CISBP-RNA/Mus_musculus.meme | cut -d" " -f3 >> ./motif_analysis/Motif_Comparison_Output/mus_final/mus_cyto_Related_RBP.temp
done
paste ./motif_analysis/Motif_Comparison_Output/mus_final/mus_cyto_Related_RBP.temp ./motif_analysis/Motif_Comparison_Output/mus_final/mus_cyto_favour_motif_tomtom.tsv > ./motif_analysis/Motif_Comparison_Output/mus_final/mus_cyto_favour_motif_final.tsv
rm ./motif_analysis/Motif_Comparison_Output/mus_final/mus_cyto_Related_RBP.temp