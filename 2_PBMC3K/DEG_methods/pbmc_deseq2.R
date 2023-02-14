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

load("counts_lognorm_scexperi.RData")
load("pbmc_label.RData")
cluster.label[c("CGGGCATGACCCAA-1","CTTGATTGATCTTC-1")] <- "Platelet"
pbmc <- pbmc[,cluster.label!="Platelet"]
cluster.label <- factor(cluster.label[cluster.label!="Platelet"])
# levels(cluster.label) <- c(1,2,1,3,1,2,1,2)

#pbmc <- as.SingleCellExperiment(pbmc)
colData(pbmc) <- cbind(colData(pbmc),cluster.label)
names(pbmc@colData@listData)[6] <- "condition"

snowparam <- SnowParam(workers = 12, type = "SOCK")
register(snowparam, default = TRUE)
registered()

colData(pbmc)[,6] <- factor(colData(pbmc)[,6])
levels(colData(pbmc)[,6]) <- 1:length(levels(colData(pbmc)[,6]))
DE.counts <- DESeqDataSetFromMatrix(pbmc@assays@data@listData$counts,colData(pbmc),~condition)
DEseq.results <- DESeq(DE.counts,parallel = T,quiet = T)
DEseq.results <- results(DEseq.results,alpha = 0.05)
# DEseq.results <- cbind(DEseq.results$padj,DEseq.results$pvalue)
save(DEseq.results,file = "./pbmc_DEG/pbmc_DEseq2.RData")