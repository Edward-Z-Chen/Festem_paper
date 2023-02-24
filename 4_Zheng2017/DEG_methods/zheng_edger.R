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

# Batch 1 -----------------------------------------------------------------

pbmc <- readRDS("./results/b1.rds")
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]
cluster.label <- factor(pbmc@meta.data$celltype)
levels(cluster.label) <- 1:nlevels(cluster.label)

pbmc <- as.SingleCellExperiment(pbmc)
names(pbmc@colData@listData)[4] <- "condition"

# snowparam <- SnowParam(workers = 12, type = "SOCK")
# register(snowparam, default = TRUE)
# registered()

colData(pbmc)[,4] <- factor(colData(pbmc)[,4])
levels(colData(pbmc)[,4]) <- 1:length(levels(colData(pbmc)[,4]))
exprSet <- DGEList(counts = as.matrix(pbmc@assays@data@listData$counts), group = colData(pbmc)[,4])
design <- model.matrix(~cluster.label)
colnames(design) <- levels(cluster.label)

if (!full_flag){
  keep <- filterByExpr(exprSet, design)
  exprSet <- exprSet[keep, , keep.lib.sizes=FALSE]
}

exprSet <- calcNormFactors(exprSet)
exprSet <- estimateDisp(exprSet,design)
# exprSet <- estimateCommonDisp(exprSet)
# exprSet <- estimateTagwiseDisp(exprSet)
# et <- exactTest(exprSet)
fit <- glmQLFit(exprSet, design)
qlf <- glmQLFTest(fit, coef=2:(nlevels(cluster.label)))
EdgeR.result <- qlf$table$PValue
names(EdgeR.result) <- rownames(qlf@.Data[[2]])
if (full_flag){
  save(EdgeR.result,file = "./results/b1_edgeR_full.RData")
}else {
  save(EdgeR.result,file = "./results/b1_edgeR.RData")
}

# Batch 2 -----------------------------------------------------------------

pbmc <- readRDS("./results/b2.rds")
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]
cluster.label <- factor(pbmc@meta.data$celltype)
levels(cluster.label) <- 1:nlevels(cluster.label)

pbmc <- as.SingleCellExperiment(pbmc)
names(pbmc@colData@listData)[4] <- "condition"

# snowparam <- SnowParam(workers = 12, type = "SOCK")
# register(snowparam, default = TRUE)
# registered()

colData(pbmc)[,4] <- factor(colData(pbmc)[,4])
levels(colData(pbmc)[,4]) <- 1:length(levels(colData(pbmc)[,4]))
exprSet <- DGEList(counts = as.matrix(pbmc@assays@data@listData$counts), group = colData(pbmc)[,4])
design <- model.matrix(~cluster.label)
colnames(design) <- levels(cluster.label)

if (!full_flag){
  keep <- filterByExpr(exprSet, design)
  exprSet <- exprSet[keep, , keep.lib.sizes=FALSE]
}

exprSet <- calcNormFactors(exprSet)
exprSet <- estimateDisp(exprSet,design)
# exprSet <- estimateCommonDisp(exprSet)
# exprSet <- estimateTagwiseDisp(exprSet)
# et <- exactTest(exprSet)
fit <- glmQLFit(exprSet, design)
qlf <- glmQLFTest(fit, coef=2:(nlevels(cluster.label)))
EdgeR.result <- qlf$table$PValue
names(EdgeR.result) <- rownames(qlf@.Data[[2]])
if (full_flag){
  save(EdgeR.result,file = "./results/b2_edgeR_full.RData")
}else {
  save(EdgeR.result,file = "./results/b2_edgeR.RData")
}