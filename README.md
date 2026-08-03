# circExor

#### Introduction
Certain circular RNAs (circRNAs) are selectively enriched in extracellular vesicles (EVs), where they contribute to intercellular communication and represent promising biomarkers, yet the sequence determinants of their sorting remain unclear. Existing computational predictors are optimized mainly for linear RNAs and rarely address circRNA localization into EVs. Here we introduce circExor, the first framework specifically designed for circRNA EV localization. We curated a dedicated benchmark dataset of 2,102 circRNAs and implemented a variable-length end-to-end concatenation strategy together with k-mer frequency encoding to accommodate circular topology, long sequence length, and length heterogeneity. Using a Tress based classifier, circExor achieved superior performance compared with RNAlocate-v3 and exoGRU, reaching an AUROC of 0.743 and AUPRC of 0.784 on the held-out test set. SHAP-based analysis, sequence perturbation experiment, and motif mapping revealed that a small number of short high-impact motifs drive classification, and implicated RBPs such as YBX1, hnRNPK, HNRNPL, and NOVA2 in circRNA sorting. circExor thus provides not only a predictive tool but also an interpretable framework that links in silico modeling to mechanistic hypotheses, supporting biomarker discovery and therapeutic design.

![Graphic Abstract](./_plot/graphic%20abstract.svg)


#### Repo Structure
This repository provides all training and evaluation code for the circExoer project, along with some training materials.

1.  `/SHAP/`: Stores scripts and results for SHAPley value calculation and visualization of machine learning results
2.  `/_plot/`: Contains images used in markdown documentation
3.  `/models/`: Code for building and training machine learning models&deep learning models, as well as saved machine learning results and outputs. Model training related analysis such as Ablation analysis stored here as well
4.  `/resources/`: Initially screened circRNA localization information from the locate database
5.  `/sample_preprocessing/`: Scripts for preprocessing resources in `/resources/` and human genome fasta files
6.  `/benchmark_evaluation/`: Comparison of circExor with other published models on baseline
7.  `/motif_analysis/`: Scripts and results related to RBP finding and motif enrichment analysis

#### Data Access

The datasets generated and processed during this study are available in this repository and in the associated Zenodo archive (https://doi.org/10.5281/zenodo.17448292). The main data files are organized as follows.

| Data item | Repository path | Description |
| --- | --- | --- |
| Curated circRNA benchmark dataset | [`sample_preprocessing/circRNA/output_with_sequences.csv`](./sample_preprocessing/circRNA/output_with_sequences.csv) | Final curated circRNA localization table used for model development. The file contains 2,102 circRNAs with RNALocate metadata, binary labels (`tag`; 1 = EV-associated, 0 = cellular-associated), circRNA sequences, and sequence lengths. |
| EV-associated benchmark sequences | [`sample_preprocessing/circRNA/EV_sequences.fasta`](./sample_preprocessing/circRNA/EV_sequences.fasta) | FASTA sequences for the 1,187 EV-associated circRNAs in the curated benchmark dataset. |
| Cellular-associated benchmark sequences | [`sample_preprocessing/circRNA/Cyto_sequences.fasta`](./sample_preprocessing/circRNA/Cyto_sequences.fasta) | FASTA sequences for the 915 cellular-associated circRNAs in the curated benchmark dataset. |
| Benchmark preprocessing workflow | [`sample_preprocessing/circRNA/circRNA_location_sequence_concentrationg.ipynb`](./sample_preprocessing/circRNA/circRNA_location_sequence_concentrationg.ipynb) | Notebook used to curate, filter, annotate, and export the benchmark circRNA dataset and FASTA files. |
| Model-input k-mer feature generation | [`models/ML_models/circRNA_EV_cell_classification_tridivided_intra5fold.py`](./models/ML_models/circRNA_EV_cell_classification_tridivided_intra5fold.py) | Script that generates the 3-mer, 4-mer, and 5-mer count and normalized-frequency feature matrices from the curated benchmark dataset before model training and validation. |
| Model-input k-mer frequency matrix | [`derived_features/circRNA_benchmark_kmer345_frequency_matrix.tsv.gz`](./derived_features/circRNA_benchmark_kmer345_frequency_matrix.tsv.gz) | Derived 2,102 x 1,344 normalized-frequency matrix for all 3-mer, 4-mer, and 5-mer features generated from the curated benchmark dataset. |
| Model-input k-mer raw-count matrix | [`derived_features/circRNA_benchmark_kmer345_rawcount_matrix.tsv.gz`](./derived_features/circRNA_benchmark_kmer345_rawcount_matrix.tsv.gz) | Derived 2,102 x 1,344 raw-count matrix for the same 3-mer, 4-mer, and 5-mer feature set. |
| SHAP value matrix | [`SHAP/circRNA_shap_output/circRNA_shap_df.txt`](./SHAP/circRNA_shap_output/circRNA_shap_df.txt) | Derived SHAP value matrix for 1,098 training-set circRNAs across 1,344 k-mer features. |
| Mean absolute SHAP scores | [`SHAP/circRNA_shap_output/circRNA_mean_abs_shap.txt`](./SHAP/circRNA_shap_output/circRNA_mean_abs_shap.txt) | Feature-level mean absolute SHAP scores for the 1,344 k-mer features. |
| Quantified k-mer importance table | [`SHAP/circRNA_shap_output/kmer_importance_quantification_info.txt`](./SHAP/circRNA_shap_output/kmer_importance_quantification_info.txt) | Summary table linking k-mer importance, SHAP direction, and localization-associated feature interpretation. |
| SHAP analysis workflow | [`SHAP/circRNA_shap.ipynb`](./SHAP/circRNA_shap.ipynb) | Notebook used to calculate and visualize SHAP-based feature importance. |
| Sequence perturbation FASTA files | [`SHAP/test_EV_sequences_mutated.fasta`](./SHAP/test_EV_sequences_mutated.fasta), [`SHAP/test_Cyto_sequences_mutated.fasta`](./SHAP/test_Cyto_sequences_mutated.fasta) | Mutated test-set circRNA sequences used for sequence perturbation analysis. |
| External breast cancer cell-line circRNA table | [`sample_preprocessing/circRNA_external_testset/BRCA_celline_circRNA_arrary.csv`](./sample_preprocessing/circRNA_external_testset/BRCA_celline_circRNA_arrary.csv) | Processed external circRNA abundance table used for external validation and sequence matching. |
| Benchmark evaluation outputs | [`benchmark_evaluation/`](./benchmark_evaluation/) | Input FASTA files, prediction files, and notebooks used to compare circExor with exoGRU and RNAlocate-v3. |
| Motif analysis outputs | [`motif_analysis/`](./motif_analysis/) | Consensus k-mer, assembled sequence, motif comparison, and enrichment-analysis files used for RBP motif interpretation. |

Large downloaded reference resources and temporary intermediate files are not all stored directly in the Git repository. In particular, large circBank/exoRBase downloads and MMseqs2 temporary sequence-comparison outputs are excluded from version control and are either regenerated by the preprocessing scripts or provided through the Zenodo archive.

#### Environment Setup

This project's dependencies are managed with Conda.

1.  **Prerequisite**: Ensure [Conda](https://docs.conda.io/en/latest/miniconda.html) is installed.

2.  **Create and activate the environment** using the `spec.txt` file:
    ````bash
    # Create the environment
    conda create --name your_env_name --file spec.txt
    
    # Activate the environment
    conda activate your_env_name
    ````

3.  **For Jupyter Notebooks**: Since most scripts are Jupyter notebooks, you may need to install `ipykernel` to make the environment available as a Jupyter kernel.
    ````bash
    conda install -n your_env_name ipykernel
    ````
#### To Do List
- [x] Gradually modify absolute paths to relative paths
- [x] Update English version annotations
- Snakemake or directly callable packaged interface

## License

This project is licensed under the MIT License - see the [LICENSE](./LICENSE) file for details.

This repository contains modified code from the [RNAlight](https://github.com/YangLab/RNAlight.git) project, which is licensed under the MIT License.
