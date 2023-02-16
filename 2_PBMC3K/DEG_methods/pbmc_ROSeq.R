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
library(ROSeq)
library(limma)
set.seed(321)
pbmc <- readRDS("./results/pbmc3k.rds")
load("./results/pbmc_label.RData")
cluster.label[c("CGGGCATGACCCAA-1","CTTGATTGATCTTC-1")] <- "Platelet"
pbmc <- pbmc[,cluster.label!="Platelet"]
cluster.label <- factor(cluster.label[cluster.label!="Platelet"])
levels(cluster.label) <- 1:nlevels(cluster.label)
counts2<-limma::voom(ROSeq::TMMnormalization(as.matrix(pbmc@assays$RNA@counts)))

cl <- makeCluster(getOption("cl.cores", 12))
roseq.tmp <- matrix(nrow = nlevels(cluster.label),ncol = nrow(pbmc))
for (l in 1:nlevels(cluster.label)){
  condition.tmp <- cluster.label == l
  condition.tmp <- as.numeric(condition.tmp)+1
  # Only 6 cores are used due to memory limits
  roseq.tmp[l,] <- ROSeq(countData=counts2$E, condition = condition.tmp, numCores=6)[,"pVals"]
}
roseq.tmp <- apply(roseq.tmp,2,min,na.rm = T)
roseq.tmp <- roseq.tmp * nlevels(cluster.label)
save(roseq.tmp,file = "./results/pbmc_ROSeq.RData")
stopCluster(cl)