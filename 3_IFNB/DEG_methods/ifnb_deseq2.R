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

ifnb <- readRDS("./results/ifnb_ctrl.rds")
cluster.labels <- ifnb@meta.data$seurat_annotations
ifnb <- as.SingleCellExperiment(ifnb)
colData(ifnb) <- cbind(colData(ifnb),cluster.labels)
names(ifnb@colData@listData)[ncol(colData(ifnb))] <- "condition"

snowparam <- SnowParam(workers = 12, type = "SOCK")
register(snowparam, default = TRUE)
registered()

colData(ifnb)[,ncol(colData(ifnb))] <- factor(colData(ifnb)[,ncol(colData(ifnb))])
levels(colData(ifnb)[,ncol(colData(ifnb))]) <- 1:length(levels(colData(ifnb)[,ncol(colData(ifnb))]))
DE.counts <- DESeqDataSetFromMatrix(ifnb@assays@data@listData$counts,colData(ifnb),~condition)
DEseq.results <- DESeq(DE.counts,parallel = T,quiet = T)
DEseq.results <- results(DEseq.results,alpha = 0.05)
# DEseq.results <- cbind(DEseq.results$padj,DEseq.results$pvalue)
save(DEseq.results,file = "./results/ifnb_DEseq2.RData")