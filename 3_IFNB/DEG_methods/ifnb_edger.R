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

ifnb <- readRDS("./results/ifnb_ctrl.rds")
cluster.labels <- ifnb@meta.data$seurat_annotations
ifnb <- as.SingleCellExperiment(ifnb)
colData(ifnb) <- cbind(colData(ifnb),cluster.labels)
names(ifnb@colData@listData)[ncol(colData(ifnb))] <- "condition"

colData(ifnb)[,ncol(colData(ifnb))] <- factor(colData(ifnb)[,ncol(colData(ifnb))])
levels(colData(ifnb)[,ncol(colData(ifnb))]) <- 1:length(levels(colData(ifnb)[,ncol(colData(ifnb))]))
exprSet <- DGEList(counts = as.matrix(ifnb@assays@data@listData$counts), group = colData(ifnb)[,ncol(colData(ifnb))])
design <- model.matrix(~cluster.labels)
colnames(design) <- levels(cluster.labels)

if (!full_flag){
  keep <- filterByExpr(exprSet, design)
  exprSet <- exprSet[keep, , keep.lib.sizes=FALSE]
}

exprSet <- calcNormFactors(exprSet)
exprSet <- estimateDisp(exprSet,design)
fit <- glmQLFit(exprSet, design)
qlf <- glmQLFTest(fit, coef=2:(nlevels(cluster.labels)))
EdgeR.result <- qlf$table$PValue
names(EdgeR.result) <- rownames(qlf@.Data[[2]])
if (full_flag){
  save(EdgeR.result,file = "./results/ifnb_edgeR_full.RData")
}else {
  save(EdgeR.result,file = "./results/ifnb_edgeR.RData")
}