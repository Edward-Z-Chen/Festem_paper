library(Seurat)
ifnb <- readRDS("./results/ifnb_ctrl.rds")
ifnb <- NormalizeData(ifnb)


# Festem ------------------------------------------------------------------
load("./results/ifnb_ctrl_Festem.RData")
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

# HVGvst ------------------------------------------------------------------
ifnb <- FindVariableFeatures(ifnb, selection.method = "vst", nfeatures = 8000,verbose = FALSE)
hvgvst <- VariableFeatures(ifnb)

# HVGdisp ------------------------------------------------------------------
ifnb <- FindVariableFeatures(ifnb, selection.method = "disp", nfeatures = 8000,verbose = FALSE)
hvgdisp <- VariableFeatures(ifnb)


# DUBStepR ----------------------------------------------------------------
library(DUBStepR)
dub.list <- DUBStepR(ifnb@assays$RNA@data)
dub <- dub.list[["optimal.feature.genes"]]
rm(dub.list)


# devianceFS --------------------------------------------------------------
library(scry)
devianceFS_out <- scry::devianceFeatureSelection(object = ifnb@assays$RNA@counts)
names(devianceFS_out) <- rownames(ifnb)
devianceFS <- names(sort(devianceFS_out, decreasing = TRUE))


# TrendVar ----------------------------------------------------------------
source("../utils/trendVar.R")
library(SingleCellExperiment)
trendvar <- trendVarFS(ifnb@assays$RNA@counts,ifnb@assays$RNA@data)
sum(trendvar[["var.out"]]@listData$FDR < 0.05,na.rm = T)
trendvar <- trendvar[["genes"]]

save(EM,hvgvst,hvgdisp,dub,devianceFS,trendvar,file = "./results/ifnb_ctrl_hvggenes.RData")