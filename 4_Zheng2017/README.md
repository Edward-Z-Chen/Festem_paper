# Analysis of the IFNB dataset

- 4.1 Preprocessing and pre-clustering. Sequencing data should be first downloaded from xxx and placed in this folder.
- 4.2 Run Festem. Cost 124.8 and 128.3 seconds (wall clock time) for each batch with 12 cores.
- 4.3 Run competing DEG detection methods.
- 4.4 Run TN test.
- 4.5 Run competing feature selection methods.
- 3.6 Clustering and UMAP reduction based on top 2500 genes selected by each method and the top 25 PCs.
- 3.7 Calculate CH indices of the PCs and clustering results derived from genes selected by feature selection methods.
- 3.8 Construct silver standard of DEGs based on housekeeping genes downloaded from https://housekeeping.unicamp.br/Housekeeping_GenesHuman.RData.
- 3.9 Festem algorithm and clustering of the stimulated group with the same procedure.
- 3.10 Figure 3A (top right), Figure 4E-G, Figure S5 (second from top), Figure S6 (top right), Figure S9, Figure S10, Figure S11, Figure S12 and Figure S13. In Figure S12A and B, we performed integrated analysis of the control and stimulated groups with *Harmony*. P-values of features were combined with the Bonferroni procedure and FDR were controlled with the Benjamini-Hochberg procedure.

**Remark**: For results in 3.3 and 3.4, we only provide a tidy version named "ifnb_ctrl_DEG_results.RData" derived from 3.10 in the "results" folder.

## Reference

Zheng, G.X.Y. et al. Massively parallel digital transcriptional profiling of single cells. Nature Communications 8 (2016)

Tran, H., Ang, K.S., Chevrier, M., Zhang, X., Lee, N.Y., Goh, M., & Chen, J. A benchmark of batch-effect correction methods for single-cell RNA sequencing data. Genome Biology, 21 (2020).

## Prerequisite

- Download housekeeping gene list from https://housekeeping.unicamp.br/Housekeeping_GenesHuman.RData and place it in **this** folder.

## Dependency

Besides the packages listed in "1_Simulation", the following R packages are also needed:

SeuratData==0.2.1, Harmony==0.1.0, plyr==1.8.6, pROC==1.17.0.1, ggsci==2.9, cowplot==1.1.1, ggalluvial==0.12.3