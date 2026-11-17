# Software setup

The root `requirements.txt` covers the numerical calculations and plots. `processing_requirements.txt` adds Python packages imported by the assay-processing code. `R_packages.txt` lists the R/Bioconductor packages. Recorded ChIP tool versions are in `lamin_b1/chipseq/software.tsv`.

R environments used for the analyses:

| Analysis | R | Packages |
|---|---|---|
| Regional B1 models | 4.1.2 | limma 3.50.3 |
| B1 and LAD-associated gene ontology | 4.2.3 | clusterProfiler 4.6.2, org.Mm.eg.db 3.16.0, GO.db 3.16.0 |
| TE differential expression and effect shrinkage | 4.3.3 | DESeq2 1.42.0, apeglm 1.24.0 |

Set the Python, R and Bash commands in root `config.toml`. Use `runtime.rscript` in a separate configuration supplied through `--config` to select the R environment for each analysis. `python run.py check` reports prerequisites before execution. The launcher does not install software.

Raw processing uses the tools named by each dataset guide, including Bismark/HISAT2, Trim Galore/Cutadapt, FastQC/MultiQC, Bowtie2, samtools, bedtools, deepTools, Salmon and Telescope 1.0.3.

Shared-HMM training uses pomegranate 0.11.1 in its original Python environment. Select that executable separately; a current pomegranate installation is not a replacement. The fitted-model decoding target does not require pomegranate.

