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
ifnb <- readRDS("./results/ifnb_ctrl.rds")
ifnb <- NormalizeData(ifnb)
cluster.label <- factor(ifnb@meta.data$seurat_annotations)
counts <- ifnb@assays$RNA@data
rm(ifnb)
counts <- as.matrix(counts)

levels(cluster.label) <- 1:nlevels(cluster.label)
cluster.labels <- cluster.label
wil.result <- vector("list",nlevels(cluster.labels))
my.wilcox <- function(x,label){
  # labels should be T or F
  wilcox.test(x[label],x[!label])$p.value
}
cl <- makeCluster(getOption("cl.cores", 12))
for (i in 1:nlevels(cluster.labels)){
  wil.result[[i]] <- parApply(cl,counts,1,my.wilcox,label = (cluster.labels==i))
}
save(wil.result,file = "./results/ifnb_wilcoxon.RData")