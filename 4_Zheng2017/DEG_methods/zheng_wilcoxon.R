library(SingleCellExperiment)
library(scater)
library(stringr)
library(stats)
library(BiocParallel)
library(BiocGenerics)
library(parallel)
library(MASS)
library(nloptr)
library(edgeR)
library(biclust)
library(Seurat)
# Batch 1 -----------------------------------------------------------------

pbmc <- readRDS("./results/b1.rds")
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]
cluster.label <- factor(pbmc@meta.data$celltype)
levels(cluster.label) <- 1:nlevels(cluster.label)
pbmc <- NormalizeData(pbmc)
counts <- pbmc@assays$RNA@data
rm(pbmc)
counts <- as.matrix(counts)

wil.result <- vector("list",nlevels(cluster.label))
my.wilcox <- function(x,label){
  # labels should be T or F
  wilcox.test(x[label],x[!label])$p.value
}
cl <- makeCluster(getOption("cl.cores", 12))
for (i in 1:nlevels(cluster.label)){
  wil.result[[i]] <- parApply(cl,counts,1,my.wilcox,label = (cluster.label==i))
}
save(wil.result,file = "./results/b1_wilcoxon.RData")

# Batch 2 -----------------------------------------------------------------

pbmc <- readRDS("./results/b2.rds")
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]
cluster.label <- factor(pbmc@meta.data$celltype)
levels(cluster.label) <- 1:nlevels(cluster.label)
pbmc <- NormalizeData(pbmc)
counts <- pbmc@assays$RNA@data
rm(pbmc)
counts <- as.matrix(counts)

wil.result <- vector("list",nlevels(cluster.label))
my.wilcox <- function(x,label){
  # labels should be T or F
  wilcox.test(x[label],x[!label])$p.value
}
cl <- makeCluster(getOption("cl.cores", 12))
for (i in 1:nlevels(cluster.label)){
  wil.result[[i]] <- parApply(cl,counts,1,my.wilcox,label = (cluster.label==i))
}
save(wil.result,file = "./results/b2_wilcoxon.RData")