library(SingleCellExperiment)
library(scater)
library(stringr)
library(stats)
library(BiocParallel)
library(BiocGenerics)
library(parallel)
library(MASS)
library(nloptr)
library(Seurat)
library(edgeR)
library(singleCellHaystack)
set.seed(321)
# Batch 1 -----------------------------------------------------------------

pbmc <- readRDS("./results/b1.rds")
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:30)
haystack_umap <- haystack(pbmc,coord = "umap",method = "2D")
haystack_pca <- haystack(pbmc,coord = "pca",dims = 1:30,method = "highD")
save(haystack_pca,haystack_umap,file = "./results/b1_haystack.RData")

pbmc <- RunUMAP(pbmc, dims = 1:20)
haystack_umap <- haystack(pbmc,coord = "umap",method = "2D")
haystack_pca <- haystack(pbmc,coord = "pca",dims = 1:20,method = "highD")
save(haystack_pca,haystack_umap,file = "./results/b1_haystack_20PC.RData")

# Batch 2 -----------------------------------------------------------------

pbmc <- readRDS("./results/b2.rds")
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:30)
haystack_umap <- haystack(pbmc,coord = "umap",method = "2D")
haystack_pca <- haystack(pbmc,coord = "pca",dims = 1:30,method = "highD")
save(haystack_pca,haystack_umap,file = "./results/b2_haystack.RData")

pbmc <- RunUMAP(pbmc, dims = 1:20)
haystack_umap <- haystack(pbmc,coord = "umap",method = "2D")
haystack_pca <- haystack(pbmc,coord = "pca",dims = 1:20,method = "highD")
save(haystack_pca,haystack_umap,file = "./results/b2_haystack_20PC.RData")