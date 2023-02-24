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

# Batch 1 -----------------------------------------------------------------

pbmc <- readRDS("./results/b1.rds")
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]
cluster.label <- factor(pbmc@meta.data$celltype)
levels(cluster.label) <- 1:nlevels(cluster.label)

pbmc <- as.SingleCellExperiment(pbmc)
names(pbmc@colData@listData)[4] <- "condition"
pbmc <- logNormCounts(pbmc)
assayNames(pbmc)[2] <- "normcounts"

snowparam <- SnowParam(workers = 12, type = "SOCK")
register(snowparam, default = TRUE)
registered()

colData(pbmc)[,4] <- factor(colData(pbmc)[,4])
levels(colData(pbmc)[,4]) <- 1:length(levels(colData(pbmc)[,4]))

mast.pbmc <- FromMatrix(as.matrix(normcounts(pbmc)),data.frame(condition=colData(pbmc)[,4],wellKey = colnames(normcounts(pbmc))),check_sanity = F)
zlm.output <- zlm(~ condition,mast.pbmc)
zlm.lr <- lrTest(zlm.output, 'condition')
mast.results <- zlm.lr[,,'Pr(>Chisq)']
rm(mast.pbmc,zlm.output,zlm.lr)

save(mast.results,file = "./results/b1_mast.RData")