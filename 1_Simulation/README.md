# Simulation

- 1.1 Estimate the parameters of the negative binomial distribution from the PBMC3K dataset, which can be downloaded from https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz. Labels are generated according to https://satijalab.org/seurat/articles/pbmc3k_tutorial.html. (pbmc_tutorial.rds upload??)

- 1.2 & 1.3 Five cell types setting with high noise (1.2) and low noise (1.3). Results are stored in the "results" folder.

## Dependency
The scripts in this repository are excuted in R (v4.0.5) and use the following R packages:

EMDE (see below), Seurat==4.0.1, BiocParallel==1.24.1, BiocGenerics==0.36.1, ggplot2==3.3.5, ggpubr==0.4.0, pracma==2.4.2, lcmix==0.3, Matrix==1.4-0, SingleCellExperiment==1.12.0, scater==1.18.6, stringr==1.4.0, MASS==7.3-55, nloptr==1.2.2.2, tidyverse==1.3.1, peakRAM==1.0.2

Competing DEG detection methods use the following R packages:

  - Wilcoxon rank sum test: wilcox.test in stat==4.0.5
  - DESeq2: DESeq2==1.30.1
  - EdgeR: edgeR==3.32.1
  - MAST: MAST==1.16.0
  - singleCellHaystack: singleCellHaystack==0.3.4
  - ROSeq: ROSeq==1.2.10, limma==3.46.0

Most of the above can be installed from CRAN and bioconductor. 

For EMDE (a preview version of Festem), install from !!!(.tar.gz). 

For lcmix, install from R-Forge:
```
install.packages("lcmix", repos="http://R-Forge.R-project.org")
```


