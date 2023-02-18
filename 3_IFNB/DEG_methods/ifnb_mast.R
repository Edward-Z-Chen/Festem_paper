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
library(MAST)

ifnb <- readRDS("./results/ifnb_ctrl.rds")
cluster.labels <- ifnb@meta.data$seurat_annotations
ifnb <- as.SingleCellExperiment(ifnb)
colData(ifnb) <- cbind(colData(ifnb),cluster.labels)
ifnb <- logNormCounts(ifnb)
assayNames(ifnb)[2] <- "normcounts"
names(ifnb@colData@listData)[ncol(colData(ifnb))] <- "condition"

snowparam <- SnowParam(workers = 12, type = "SOCK")
register(snowparam, default = TRUE)
registered()

colData(ifnb)[,ncol(colData(ifnb))] <- factor(colData(ifnb)[,ncol(colData(ifnb))])
levels(colData(ifnb)[,ncol(colData(ifnb))]) <- 1:length(levels(colData(ifnb)[,ncol(colData(ifnb))]))

mast.ifnb <- FromMatrix(as.matrix(normcounts(ifnb)),data.frame(condition=colData(ifnb)[,ncol(colData(ifnb))],wellKey = colnames(normcounts(ifnb))),check_sanity = F)
zlm.output <- zlm(~ condition,mast.ifnb)
zlm.lr <- lrTest(zlm.output, 'condition')
mast.results <- zlm.lr[,,'Pr(>Chisq)']
rm(mast.ifnb,zlm.output,zlm.lr)

save(mast.results,file = "./results/ifnb_mast.RData")