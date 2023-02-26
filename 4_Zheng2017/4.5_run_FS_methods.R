library(SingleCellExperiment)
library(scater)
library(stringr)
library(Seurat)
library(stats)
library(BiocParallel)
library(BiocGenerics)
library(parallel)
library(MASS)
library(nloptr)
library(edgeR)
library(ggplot2)
library(R.utils)
# Batch 1 -----------------------------------------------------------------

pbmc <- readRDS("./results/b1.rds")
pbmc <- NormalizeData(pbmc)
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]

load("./results/b1_Festem.RData")
## EM
result.EM = data.frame(names = colnames(em.result), 
                       p = p.adjust(em.result[1,],"BH"), 
                       EM = em.result.9[2,])
tmp.na <- result.EM[is.na(result.EM$p),1]
result.EM <- result.EM[!is.na(result.EM$p),]
gene.names <- result.EM[result.EM$p<0.05 & result.EM$EM>0,]
gene.names <- gene.names[order(-gene.names$EM),]
tmp <- result.EM[result.EM$p>=0.05 & result.EM$EM>0,]
tmp <- tmp[order(-tmp$EM),]
gene.names <- rbind(gene.names,tmp)
tmp <- result.EM[result.EM$EM<=0,]
tmp <- tmp[order(tmp$p,-tmp$EM),]
gene.names <- rbind(gene.names,tmp)
gene.names <- gene.names[,1]
gene.names <- c(gene.names,tmp.na)
EM <- gene.names

## HVGvst
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 8000,verbose = FALSE)
hvgvst <- VariableFeatures(pbmc)

## HVGDisp
tryCatch({
  withTimeout({
    pbmc <- FindVariableFeatures(object = pbmc, selection.method = "disp", nfeatures = 8000,verbose = FALSE)
    hvgdisp <- VariableFeatures(pbmc)
  },timeout = 30000,onTimeout = "error")
}, error = function(e){
  hvgdisp <- NULL
  cat("HVGdisp Error!\n")
  cat("ERROR :",conditionMessage(e),"\n")
})

## DUBStepR
library(DUBStepR)
dub.list <- DUBStepR(pbmc@assays$RNA@data)
dub <- dub.list[["optimal.feature.genes"]]
rm(dub.list)


##devianceFS
library(scry)
devianceFS_out <- scry::devianceFeatureSelection(object = pbmc@assays$RNA@counts)
names(devianceFS_out) <- rownames(pbmc)
devianceFS <- names(sort(devianceFS_out, decreasing = TRUE))
rm(devianceFS_out)

## TrendVar
source("../utils/trendVar.R")
library(SingleCellExperiment)
trendvar <- trendVarFS(pbmc@assays$RNA@counts,pbmc@assays$RNA@data)
sum(trendvar[["var.out"]]@listData$FDR < 0.05,na.rm = T)
trendvar <- trendvar[["genes"]]
 
save(EM,hvgvst,hvgdisp,dub,devianceFS,trendvar,file = "./results/b1_hvggenes.RData")

# Batch 2 -----------------------------------------------------------------

pbmc <- readRDS("./results/b2.rds")
pbmc <- NormalizeData(pbmc)
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]
load("./results/b2_Festem.RData")
## EM
result.EM = data.frame(names = colnames(em.result), 
                       p = p.adjust(em.result[1,],"BH"), 
                       EM = em.result.9[2,])
tmp.na <- result.EM[is.na(result.EM$p),1]
result.EM <- result.EM[!is.na(result.EM$p),]
gene.names <- result.EM[result.EM$p<0.05 & result.EM$EM>0,]
gene.names <- gene.names[order(-gene.names$EM),]
tmp <- result.EM[result.EM$p>=0.05 & result.EM$EM>0,]
tmp <- tmp[order(-tmp$EM),]
gene.names <- rbind(gene.names,tmp)
tmp <- result.EM[result.EM$EM<=0,]
tmp <- tmp[order(tmp$p,-tmp$EM),]
gene.names <- rbind(gene.names,tmp)
gene.names <- gene.names[,1]
gene.names <- c(gene.names,tmp.na)
EM <- gene.names

## HVGvst
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 8000,verbose = FALSE)
hvgvst <- VariableFeatures(pbmc)

## HVGDisp
tryCatch({
  withTimeout({
    pbmc <- FindVariableFeatures(object = pbmc, selection.method = "disp", nfeatures = 8000,verbose = FALSE)
    hvgdisp <- VariableFeatures(pbmc)
  },timeout = 30000,onTimeout = "error")
}, error = function(e){
  hvgdisp <- NULL
  cat("HVGdisp Error!\n")
  cat("ERROR :",conditionMessage(e),"\n")
})

## DUBStepR
library(DUBStepR)
dub.list <- DUBStepR(pbmc@assays$RNA@data)
dub <- dub.list[["optimal.feature.genes"]]
rm(dub.list)


##devianceFS
library(scry)
devianceFS_out <- scry::devianceFeatureSelection(object = pbmc@assays$RNA@counts)
names(devianceFS_out) <- rownames(pbmc)
devianceFS <- names(sort(devianceFS_out, decreasing = TRUE))
rm(devianceFS_out)

## TrendVar
source("../utils/trendVar.R")
library(SingleCellExperiment)
trendvar <- trendVarFS(pbmc@assays$RNA@counts,pbmc@assays$RNA@data)
sum(trendvar[["var.out"]]@listData$FDR < 0.05,na.rm = T)
trendvar <- trendvar[["genes"]]

save(EM,hvgvst,hvgdisp,dub,devianceFS,trendvar,file = "./results/b2_hvggenes.RData")
