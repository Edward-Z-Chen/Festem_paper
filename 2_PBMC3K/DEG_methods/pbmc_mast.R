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
pbmc <- readRDS("pbmc3k.rds")
pbmc <- as.SingleCellExperiment(pbmc)
load("./results/pbmc3k_label.RData")
cluster.label[c("CGGGCATGACCCAA-1","CTTGATTGATCTTC-1")] <- "Platelet"
pbmc <- pbmc[,cluster.label!="Platelet"]
cluster.label <- factor(cluster.label[cluster.label!="Platelet"])

colData(pbmc) <- cbind(colData(pbmc),cluster.label)
pbmc <- logNormCounts(pbmc)
assayNames(pbmc)[2] <- "normcounts"
names(pbmc@colData@listData)[ncol(colData(pbmc))] <- "condition"

snowparam <- SnowParam(workers = 12, type = "SOCK")
register(snowparam, default = TRUE)
registered()

colData(pbmc)[,ncol(colData(pbmc))] <- factor(colData(pbmc)[,ncol(colData(pbmc))])
levels(colData(pbmc)[,ncol(colData(pbmc))]) <- 1:length(levels(colData(pbmc)[,ncol(colData(pbmc))]))

mast.pbmc <- FromMatrix(as.matrix(normcounts(pbmc)),data.frame(condition=colData(pbmc)[,ncol(colData(pbmc))],wellKey = colnames(normcounts(pbmc))),check_sanity = F)
zlm.output <- zlm(~ condition,mast.pbmc)
zlm.lr <- lrTest(zlm.output, 'condition')
mast.results <- zlm.lr[,,'Pr(>Chisq)']
rm(mast.pbmc,zlm.output,zlm.lr)

save(mast.results,file = "./results/pbmc_mast.RData")