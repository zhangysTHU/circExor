#!/bin/bash
# 根据kmer和circRNA序列进行kmer组装
python kmer_assmble_MSA.py --kmers /lulabdata3/huangkeyun/zhangys/RNA_locator/archived/circExor_archived_2025_5/SHAP/cyto_kmer.fasta --circ /lulabdata3/huangkeyun/zhangys/RNA_locator/archived/circExor_archived_2025_5/sample_preprocessing/circRNA/cyto_sequence.fasta --out ./assmbled_kmers/cyto_nt.fasta --dedup
python kmer_assmble_MSA.py --kmers /lulabdata3/huangkeyun/zhangys/RNA_locator/archived/circExor_archived_2025_5/SHAP/EV_kmer.fasta --circ /lulabdata3/huangkeyun/zhangys/RNA_locator/archived/circExor_archived_2025_5/sample_preprocessing/circRNA/EV_sequence.fasta --out ./assmbled_kmers/EV_nt.fasta --dedup

# 对组装后的kmer进行多序列比对并生成共识序列
python MSA_consensus.py -i ./assmbled_kmers/cyto_nt.fasta -o ./consensus_kmers/cyto
python MSA_consensus.py -i ./assmbled_kmers/EV_nt.fasta -o ./consensus_kmers/EV

# 将公示序列通过TomTom比对到已知RBP motif
bash /lulabdata3/huangkeyun/zhangys/RNA_locator/motif_analysis/cyto_motif_comparision.sh
bash /lulabdata3/huangkeyun/zhangys/RNA_locator/motif_analysis/EV_motif_comparision.sh