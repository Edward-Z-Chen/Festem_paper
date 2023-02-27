library(Seurat)
library(DUBStepR)
load("./results/iCCA_Festem.RData")
nonepi <- readRDS("NonEpi.rds")

all.genes <- rownames(nonepi)

nonepi <- FindVariableFeatures(nonepi, selection.method = "vst", nfeatures = 10000,verbose = FALSE)
hvgvst <- VariableFeatures(nonepi)

## EM
p.tmp <- matrix(NA,nrow = length(all.genes),ncol = length(em.result))
rownames(p.tmp) <- all.genes
for (i in 1:length(em.result)){
  p.tmp[colnames(em.result[[i]]),i] <- em.result[[i]][1,]
}
my.min <- function(x){
  if (sum(is.na(x))==length(x)){
    NA
  } else {
    min(x,na.rm = T)
  }
}
p.tmp <- apply(p.tmp,1,my.min)
p.tmp <- p.tmp*length(em.result)
em.tmp <- matrix(NA,nrow = length(all.genes),ncol = length(em.result))
rownames(em.tmp) <- all.genes
for (i in 1:length(em.result)){
  em.tmp[colnames(em.result.9[[i]]),i] <- em.result.9[[i]][2,]
}
my.max <- function(x){
  if (sum(is.na(x))==length(x)){
    NA
  } else {
    max(x,na.rm = T)
  }
}
em.tmp <- apply(em.tmp,1,my.max)

result.EM = data.frame(names = all.genes, p = p.adjust(p.tmp,"BH"), EM = em.tmp)
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

## HVGDisp10000
nonepi <- FindVariableFeatures(object = nonepi, selection.method = "disp", nfeatures = 10000,verbose = FALSE)
hvgdisp <- VariableFeatures(nonepi)

## DUBStepR
dub.list <- DUBStepR(nonepi@assays$RNA@data)
dub <- dub.list[["optimal.feature.genes"]]
rm(dub.list)

##devianceFS
library(scry)
devianceFS_out <- scry::devianceFeatureSelection(object = nonepi@assays$RNA@counts,batch = factor(nonepi@meta.data$Patient))
names(devianceFS_out) <- rownames(nonepi)
devianceFS <- names(sort(devianceFS_out, decreasing = TRUE))

## TrendVar
source("../utils/trendVar.R")
library(SingleCellExperiment)
trendVarFS <- function(counts, data,block) {
  st <- system.time({
    sce <- SingleCellExperiment(list(counts = counts, logcounts = data))
    mgvar <- scran::modelGeneVar(x = sce)
    top.hvgs <- scran::getTopHVGs(mgvar, n = nrow(mgvar))
  })
  
  return(list("var.out" = mgvar, "genes" = top.hvgs, "st" = st))
}
trendvar <- trendVarFS(nonepi@assays$RNA@counts,nonepi@assays$RNA@data,block = factor(nonepi@meta.data$Patient))
sum(trendvar[["var.out"]]@listData$FDR < 0.05,na.rm = T)
trendvar <- trendvar["genes"]

save(hvgvst,EM,hvgdisp,dub,devianceFS,trendvar,file = "./results/iCCA_hvggenes.RData")
