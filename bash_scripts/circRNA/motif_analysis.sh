#!/bin/bash
# Assemble kmers and circRNA sequences into kmer sequences
python kmer_assmble_MSA.py --kmers /lulabdata3/huangkeyun/zhangys/RNA_locator/archived/circExor_archived_2025_5/SHAP/cyto_kmer.fasta --circ /lulabdata3/huangkeyun/zhangys/RNA_locator/archived/circExor_archived_2025_5/sample_preprocessing/circRNA/cyto_sequence.fasta --out ./assmbled_kmers/cyto_nt.fasta --dedup
python kmer_assmble_MSA.py --kmers /lulabdata3/huangkeyun/zhangys/RNA_locator/archived/circExor_archived_2025_5/SHAP/EV_kmer.fasta --circ /lulabdata3/huangkeyun/zhangys/RNA_locator/archived/circExor_archived_2025_5/sample_preprocessing/circRNA/EV_sequence.fasta --out ./assmbled_kmers/EV_nt.fasta --dedup

# Perform multiple sequence alignment on assembled kmers and generate consensus sequences
python MSA_consensus.py -i ./assmbled_kmers/cyto_nt.fasta -o ./consensus_kmers/cyto
python MSA_consensus.py -i ./assmbled_kmers/EV_nt.fasta -o ./consensus_kmers/EV

# Compare consensus sequences to known RBP motifs using TomTom
bash /lulabdata3/huangkeyun/zhangys/RNA_locator/motif_analysis/cyto_motif_comparision.sh
bash /lulabdata3/huangkeyun/zhangys/RNA_locator/motif_analysis/EV_motif_comparision.sh