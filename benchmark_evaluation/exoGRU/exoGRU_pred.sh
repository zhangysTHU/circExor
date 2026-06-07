#!/bin/bash

PYTHON_EXEC="python"
SCRIPT="circExor/repeated_works/exosomians/pred_EV_fasta.py"
BOOTSTRAP_DIR="circExor/benchmark_evaluation/exoGRU/bootstrap_results"

for i in {1..5}
do
    INPUT_FASTA="${BOOTSTRAP_DIR}/bootstrap_set_${i}.fasta"
    OUTPUT_CSV="${BOOTSTRAP_DIR}/bootstrap_set_${i}_pred.csv"
    
    echo "Running prediction for bootstrap set ${i}..."
    ${PYTHON_EXEC} ${SCRIPT} ${INPUT_FASTA} ${OUTPUT_CSV}
    echo "Done. Results saved to ${OUTPUT_CSV}"
    echo "----------------------------------------"
done