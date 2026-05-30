#!/bin/bash

REF="/lulabdata3/huangkeyun/zhangys/RNA_locator/ML_python_scripts/ML_models/circRNA_ML_Model_tridivided_extra4fold_Output/train_val_set_sequences.fasta"
TARGETS=("./outputs/target_exoRbase.fasta" "./outputs/target_celline.fasta")
IDENTITY=0.9
OUT_DIR="./outputs"
TMP_DIR="$OUT_DIR/tmp_mmseqs"

mkdir -p "$TMP_DIR"

echo "Building Reference Index..."
mmseqs createdb "$REF" "$TMP_DIR/refDB" --dbtype 2
mmseqs createindex "$TMP_DIR/refDB" "$TMP_DIR" --search-type 3

for TARGET_FASTA in "${TARGETS[@]}"; do
    NAME=$(basename "$TARGET_FASTA" .fasta)
    echo "Running MMseqs2 Search for $NAME..."

    mmseqs createdb "$TARGET_FASTA" "$TMP_DIR/${NAME}_DB" --dbtype 2
    
    mmseqs search "$TMP_DIR/${NAME}_DB" "$TMP_DIR/refDB" "$TMP_DIR/${NAME}_res" "$TMP_DIR" \
        --min-seq-id "$IDENTITY" -c 0.9 --search-type 3 --threads 16

    # 导出并去重匹配到的目标ID
    mmseqs convertalis "$TMP_DIR/${NAME}_DB" "$TMP_DIR/refDB" "$TMP_DIR/${NAME}_res" "$TMP_DIR/${NAME}_hits.txt" --format-output "query"
    sort -u "$TMP_DIR/${NAME}_hits.txt" > "$OUT_DIR/${NAME}_matched_ids.lst"
    
    echo "Saved high-sim IDs to: $OUT_DIR/${NAME}_matched_ids.lst"
done