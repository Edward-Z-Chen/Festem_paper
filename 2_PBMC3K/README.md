# Analysis of the PBMC3K dataset

- 2.1 Preprocessing according to Seurat's tutorial (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) and pre-clustering. Sequencing data should be first downloaded from https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz and placed in this folder.
- 2.2 Run Festem. Cost 37.5 seconds (wall clock time) with 12 cores.
- 2.3 Run competing DEG detection methods.
- 2.4 Run TN test.
- 2.5 Run competing feature selection methods.
- 2.6 Clustering and UMAP reduction based on top 1000 genes selected by each method and the top 15 PCs.
- 2.7 Calculate CH indices of the PCs and clustering results derived from genes selected by feature selection methods.

**Remark**: For results in 2.3 and 2.4, we only provide a tidy version derived from 2.x in the "results" folder.
