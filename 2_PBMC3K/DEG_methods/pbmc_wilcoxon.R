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
pbmc <- readRDS("pbmc3k.rds")
load("./results/pbmc3k_label.RData")
cluster.label[c("CGGGCATGACCCAA-1","CTTGATTGATCTTC-1")] <- "Platelet"
pbmc <- pbmc[,cluster.label!="Platelet"]
cluster.label <- factor(cluster.label[cluster.label!="Platelet"])
counts <- pbmc@assays$RNA@data
rm(pbmc)
counts <- as.matrix(counts)

levels(cluster.label) <- 1:8
cluster.labels <- cluster.label
wil.result <- vector("list",8)
my.wilcox <- function(x,label){
  # labels should be T or F
  wilcox.test(x[label],x[!label])$p.value
}
cl <- makeCluster(getOption("cl.cores", 12))
for (i in 1:8){
  wil.result[[i]] <- parApply(cl,counts,1,my.wilcox,label = (cluster.labels==i))
}
stopCluster(cl)
save(wil.result,file = "./results/pbmc_wilcoxon.RData")