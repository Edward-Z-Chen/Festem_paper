library(SingleCellExperiment)
library(scater)
library(Seurat)
library(stringr)
library(stats)
library(BiocParallel)
library(BiocGenerics)
library(parallel)
library(MASS)
library(nloptr)
library(DESeq2)
pbmc <- readRDS("./results/pbmc3k.rds")
pbmc <- as.SingleCellExperiment(pbmc)
load("./results/pbmc3k_label.RData")
cluster.label[c("CGGGCATGACCCAA-1","CTTGATTGATCTTC-1")] <- "Platelet"
pbmc <- pbmc[,cluster.label!="Platelet"]
cluster.label <- factor(cluster.label[cluster.label!="Platelet"])

colData(pbmc) <- cbind(colData(pbmc),cluster.label)
names(pbmc@colData@listData)[ncol(colData(pbmc))] <- "condition"

snowparam <- SnowParam(workers = 12, type = "SOCK")
register(snowparam, default = TRUE)
registered()

colData(pbmc)[,ncol(colData(pbmc))] <- factor(colData(pbmc)[,ncol(colData(pbmc))])
levels(colData(pbmc)[,ncol(colData(pbmc))]) <- 1:length(levels(colData(pbmc)[,ncol(colData(pbmc))]))
DE.counts <- DESeqDataSetFromMatrix(pbmc@assays@data@listData$counts,colData(pbmc),~condition)
DEseq.results <- DESeq(DE.counts,parallel = T,quiet = T)
DEseq.results <- results(DEseq.results,alpha = 0.05)
save(DEseq.results,file = "./results/pbmc_DEseq2.RData")