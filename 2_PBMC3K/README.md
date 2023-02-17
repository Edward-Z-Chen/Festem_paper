# Analysis of the PBMC3K dataset

- 2.1 Preprocessing according to Seurat's tutorial (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) and pre-clustering. Sequencing data should be first downloaded from https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz and placed in this folder.
- 2.2 Run Festem. Cost 37.5 seconds (wall clock time) with 12 cores.
- 2.3 Run competing DEG detection methods.
- 2.4 Run TN test.
- 2.5 Run competing feature selection methods.
- 2.6 Clustering and UMAP reduction based on top 1000 genes selected by each method and the top 15 PCs.
- 2.7 Calculate CH indices of the PCs and clustering results derived from genes selected by feature selection methods.
- 2.8 Construct silver standard of DEGs based on housekeeping genes downloaded from https://housekeeping.unicamp.br/Housekeeping_GenesHuman.RData.

**Remark**: For results in 2.3 and 2.4, we only provide a tidy version derived from 2.9 in the "results" folder.

## Prerequisite

- Download preprocessed data from xxx (recommanded) and place it under the **"results"** folder, OR download the original data from https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz.
- Download housekeeping gene list from https://housekeeping.unicamp.br/Housekeeping_GenesHuman.RData and place it in **this** folder.

## Dependency

Besides the packages listed in "1_Simulation", the following R packages are also needed:

plyr==1.8.6, pROC==1.17.0.1, ggsci==2.9, cowplot==1.1.1, ggalluvial==0.12.3
