# Analysis of the IFNB dataset

- 4.1 Preprocessing and pre-clustering. Sequencing data should be first downloaded from xxx and placed in this folder. Sequencing data was downloaded from https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking (dataset 5).
- 4.2 Run Festem. Cost 124.8 and 128.3 seconds (wall clock time) for each batch with 12 cores.
- 4.3 Run competing DEG detection methods.
- 4.4 Run TN test.
- 4.5 Run competing feature selection methods.
- 4.6 Calculate CH indices of the PCs and clustering results derived from genes selected by feature selection methods.
- 4.7 Construct silver standard of DEGs based on housekeeping genes downloaded from https://housekeeping.unicamp.br/Housekeeping_GenesHuman.RData.
- 4.8 Figure 3A (bottom), Figure S5 (second from bottom and the bottom one) and Figure S6 (bottom).

**Remark**: For results in 4.3 and 4.4, we only provide a tidy version named "ifnb_ctrl_DEG_results.RData" derived from 4.8 in the "results" folder.

## Reference

Zheng, G.X.Y. et al. Massively parallel digital transcriptional profiling of single cells. Nature Communications 8 (2016)

Tran, H., Ang, K.S., Chevrier, M., Zhang, X., Lee, N.Y., Goh, M., & Chen, J. A benchmark of batch-effect correction methods for single-cell RNA sequencing data. Genome Biology, 21 (2020).

## Prerequisite

- Download the preprocessed data from xxx

- Download housekeeping gene list from https://housekeeping.unicamp.br/Housekeeping_GenesHuman.RData and place it in **this** folder.

## Dependency

Besides the packages listed in "1_Simulation", the following R packages are also needed:

SeuratData==0.2.1, plyr==1.8.6, pROC==1.17.0.1.