#!/bin/bash
python kmer_assmble_MSA.py --kmers ./SHAP/cyto_kmer.fasta --circ ./sample_preprocessing/circRNA/cyto_sequence.fasta --out ./assmbled_kmers/cyto_nt.fasta --dedup
python kmer_assmble_MSA.py --kmers ./SHAP/EV_kmer.fasta --circ ./sample_preprocessing/circRNA/EV_sequence.fasta --out ./assmbled_kmers/EV_nt.fasta --dedup

python MSA_consensus.py -i ./assmbled_kmers/cyto_nt.fasta -o ./consensus_kmers/cyto
python MSA_consensus.py -i ./assmbled_kmers/EV_nt.fasta -o ./consensus_kmers/EV

bash ./motif_analysis/cyto_motif_comparision.sh
bash ./motif_analysis/EV_motif_comparision.sh