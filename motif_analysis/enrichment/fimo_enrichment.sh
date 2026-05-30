# 1. 用 EV-Motif 扫描 EV 序列
fimo --oc fimo_EV_motif_on_EV_seq ./meme_files/EV_specific_motifs.meme /lulabdata3/huangkeyun/zhangys/RNA_locator/motif_analysis/circRNA_ligated/EV_sequences.fasta

# 2. 用 EV-Motif 扫描 Cyto 序列
fimo --oc fimo_EV_motif_on_Cyto_seq ./meme_files/EV_specific_motifs.meme /lulabdata3/huangkeyun/zhangys/RNA_locator/motif_analysis/circRNA_ligated/Cyto_sequences.fasta

# 3. 用 Cyto-Motif 扫描 EV 序列
fimo --oc fimo_Cyto_motif_on_EV_seq ./meme_files/Cyto_specific_motifs.meme /lulabdata3/huangkeyun/zhangys/RNA_locator/motif_analysis/circRNA_ligated/EV_sequences.fasta

# 4. 用 Cyto-Motif 扫描 Cyto 序列
fimo --oc fimo_Cyto_motif_on_Cyto_seq ./meme_files/Cyto_specific_motifs.meme /lulabdata3/huangkeyun/zhangys/RNA_locator/motif_analysis/circRNA_ligated/Cyto_sequences.fasta