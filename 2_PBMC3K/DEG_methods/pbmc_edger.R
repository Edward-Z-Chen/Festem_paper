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
library(edgeR)
pbmc <- readRDS("./results/pbmc3k.rds")
pbmc <- as.SingleCellExperiment(pbmc)
load("./results/pbmc3k_label.RData")
cluster.label[c("CGGGCATGACCCAA-1","CTTGATTGATCTTC-1")] <- "Platelet"
pbmc <- pbmc[,cluster.label!="Platelet"]
cluster.label <- factor(cluster.label[cluster.label!="Platelet"])

colData(pbmc) <- cbind(colData(pbmc),cluster.label)
names(pbmc@colData@listData)[ncol(colData(pbmc))] <- "condition"
# snowparam <- SnowParam(workers = 12, type = "SOCK")
# register(snowparam, default = TRUE)
# registered()

colData(pbmc)[,ncol(colData(pbmc))] <- factor(colData(pbmc)[,ncol(colData(pbmc))])
levels(colData(pbmc)[,ncol(colData(pbmc))]) <- 1:length(levels(colData(pbmc)[,ncol(colData(pbmc))]))
exprSet <- DGEList(counts = as.matrix(pbmc@assays@data@listData$counts), group = colData(pbmc)[,ncol(colData(pbmc))])
design <- model.matrix(~cluster.label)
colnames(design) <- levels(cluster.label)

if (!full_flag){
  keep <- filterByExpr(exprSet, design)
  exprSet <- exprSet[keep, , keep.lib.sizes=FALSE]
}

exprSet <- calcNormFactors(exprSet)
exprSet <- estimateDisp(exprSet,design)
fit <- glmQLFit(exprSet, design)
qlf <- glmQLFTest(fit, coef=2:(nlevels(cluster.label)))
EdgeR.result <- qlf$table$PValue
names(EdgeR.result) <- rownames(qlf@.Data[[2]])
if (full_flag){
  save(EdgeR.result,file = "./results/pbmc_edgeR_full.RData")
}else {
  save(EdgeR.result,file = "./results/pbmc_edgeR.RData")
}