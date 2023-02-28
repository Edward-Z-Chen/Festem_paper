library(Seurat)
library(BiocParallel)
library(BiocGenerics)
library(parallel)
source("../utils/type_allocation.R")
load("./results/iCCA_clustering_UMAP.RData")
gene.names <- gene.list[[1]]
nonepi <- readRDS("NonEpi.rds")
nonepi <- NormalizeData(nonepi)
nonepi@active.ident <- factor(label.list[[1]])
names(nonepi@active.ident) <- colnames(nonepi)

cl <- makeCluster(getOption("cl.cores", 20))
nonepi.gene.allocation <- parApply(cl,nonepi@assays$RNA@data[gene.names,], 1, 
                              type_allocate, type = factor(nonepi@active.ident),
                              sig.level = 0.05)
save(nonepi.gene.allocation,file = "./results/iCCA_marker_allocation.RData")
stopCluster(cl)
