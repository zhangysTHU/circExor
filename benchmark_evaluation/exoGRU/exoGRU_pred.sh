#!/bin/bash

# 定义路径变量
PYTHON_EXEC="/lulabdata3/huangkeyun/miniconda/bin/python"
SCRIPT="/lulabdata3/huangkeyun/zhangys/RNA_locator/repeated_works/exosomians/pred_EV_fasta.py"
BOOTSTRAP_DIR="/lulabdata3/huangkeyun/zhangys/RNA_locator/benchmark_evaluation/exoGRU/bootstrap_results"

# 循环 1 到 5
for i in {1..5}
do
    INPUT_FASTA="${BOOTSTRAP_DIR}/bootstrap_set_${i}.fasta"
    OUTPUT_CSV="${BOOTSTRAP_DIR}/bootstrap_set_${i}_pred.csv"
    
    echo "Running prediction for bootstrap set ${i}..."
    ${PYTHON_EXEC} ${SCRIPT} ${INPUT_FASTA} ${OUTPUT_CSV}
    echo "Done. Results saved to ${OUTPUT_CSV}"
    echo "----------------------------------------"
done