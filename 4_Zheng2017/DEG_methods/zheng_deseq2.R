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


# Batch 1 -----------------------------------------------------------------

pbmc <- readRDS("./results/b1.rds")
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]

pbmc <- as.SingleCellExperiment(pbmc)
names(pbmc@colData@listData)[4] <- "condition"

snowparam <- SnowParam(workers = 12, type = "SOCK")
register(snowparam, default = TRUE)
registered()

colData(pbmc)[,4] <- factor(colData(pbmc)[,4])
levels(colData(pbmc)[,4]) <- 1:length(levels(colData(pbmc)[,4]))
DE.counts <- DESeqDataSetFromMatrix(pbmc@assays@data@listData$counts,colData(pbmc),~condition)
DEseq.results <- DESeq(DE.counts,parallel = T,quiet = T)
DEseq.results <- results(DEseq.results,alpha = 0.05)
# DEseq.results <- cbind(DEseq.results$padj,DEseq.results$pvalue)
save(DEseq.results,file = "./results/b1_DEseq2.RData")

# Batch 2 -----------------------------------------------------------------
pbmc <- readRDS("./results/b2.rds")
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]

pbmc <- as.SingleCellExperiment(pbmc)
names(pbmc@colData@listData)[4] <- "condition"

snowparam <- SnowParam(workers = 12, type = "SOCK")
register(snowparam, default = TRUE)
registered()

colData(pbmc)[,4] <- factor(colData(pbmc)[,4])
levels(colData(pbmc)[,4]) <- 1:length(levels(colData(pbmc)[,4]))
DE.counts <- DESeqDataSetFromMatrix(pbmc@assays@data@listData$counts,colData(pbmc),~condition)
DEseq.results <- DESeq(DE.counts,parallel = T,quiet = T)
DEseq.results <- results(DEseq.results,alpha = 0.05)
# DEseq.results <- cbind(DEseq.results$padj,DEseq.results$pvalue)
save(DEseq.results,file = "./results/b2_DEseq2.RData")
