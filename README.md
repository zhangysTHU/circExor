# circExor

#### Introduction
Certain circular RNAs (circRNAs) are selectively enriched in extracellular vesicles (EVs), where they contribute to intercellular communication and represent promising biomarkers, yet the sequence determinants of their sorting remain unclear. Existing computational predictors are optimized mainly for linear RNAs and rarely address circRNA localization into EVs. Here we introduce circExor, the first framework specifically designed for circRNA EV localization. We curated a dedicated benchmark dataset of 2,102 circRNAs and implemented a variable-length end-to-end concatenation strategy together with k-mer frequency encoding to accommodate circular topology, long sequence length, and length heterogeneity. Using a Tress based classifier, circExor achieved superior performance compared with RNAlocate-v3 and exoGRU, reaching an AUROC of 0.743 and AUPRC of 0.784 on the held-out test set. SHAP-based analysis, sequence perturbation experiment, and motif mapping revealed that a small number of short high-impact motifs drive classification, and implicated RBPs such as YBX1, hnRNPK, HNRNPL, and NOVA2 in circRNA sorting. circExor thus provides not only a predictive tool but also an interpretable framework that links in silico modeling to mechanistic hypotheses, supporting biomarker discovery and therapeutic design.

![Graphic Abstract](./_plot/graphic%20abstract.svg)


#### Repo Structure
This repository provides all training and evaluation code for the circExoer project, along with some training materials.

1.  `/SHAP/`: Stores scripts and results for SHAPley value calculation and visualization of machine learning results
2.  `/_plot/`: Contains images used in markdown documentation
3.  `/models/`: Code for building and training machine learning models, as well as saved machine learning results and outputs
4.  `/resources/`: Initially screened circRNA localization information from the locate database
5.  `/sample_preprocessings/`: Scripts for preprocessing resources in `/resources/` and human genome fasta files
6.  `/references/`: circRNA fasta files obtained from public databases
7.  `/benchmark_evaluation/`: Comparison of circExor with other published models on baseline
8.  `/motif_analysis/`: Scripts and results related to RBP finding

#### To Do List
- Gradually modify absolute paths to relative paths
- Snakemake or directly callable packaged interface

## License

This project is licensed under the MIT License - see the [LICENSE](./LICENSE) file for details.

This repository contains modified code from the [RNAlight](https://github.com/YangLab/RNAlight.git) project, which is licensed under the MIT License.