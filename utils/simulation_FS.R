library(Seurat)
library(DUBStepR)

all.genes <- rownames(seurat.SD)
# Festem (Group parameter G = g_1) --------------
FS.time.mat[i,1] <- time.mat[i,2]
result.EM = data.frame(names = all.genes, p = p.adjust(em.result[1,],"BH"), EM = em.result.9[2,])
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
FS.gene.list[[1]][[i]] <- gene.names

# Festem (Group parameter G = g_2) --------------
FS.time.mat[i,10] <- time.mat[i,15]
result.EM = data.frame(names = all.genes, p = p.adjust(em.result.3g[1,],"BH"), EM = em.result.9.3g[2,])
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
FS.gene.list[[10]][[i]] <- gene.names

# Festem (Group parameter G = g_3) --------------
FS.time.mat[i,11] <- time.mat[i,17]
result.EM = data.frame(names = all.genes, p = p.adjust(em.result.4g[1,],"BH"), EM = em.result.9.4g[2,])
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
FS.gene.list[[11]][[i]] <- gene.names

# Festem (Group parameter G = g_4) --------------
FS.time.mat[i,12] <- time.mat[i,19]
result.EM = data.frame(names = all.genes, p = p.adjust(em.result.5g[1,],"BH"), EM = em.result.9.5g[2,])
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
FS.gene.list[[12]][[i]] <- gene.names

# HVGvst --------------
time.tmp <- peakRAM(
seurat.SD <- FindVariableFeatures(seurat.SD, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
)
FS.time.mat[i,2] <- time.tmp$Elapsed_Time_sec
total.memory.usage[i,"HVGvst"] <- time.tmp$Total_RAM_Used_MiB
peak.memory.usage[i,"HVGvst"] <- time.tmp$Peak_RAM_Used_MiB

FS.gene.list[[2]][[i]] <- VariableFeatures(seurat.SD)

# HVGdisp --------------
time.tmp <- peakRAM(
seurat.SD <- FindVariableFeatures(object = seurat.SD, selection.method = "disp", nfeatures = 2000,verbose = FALSE)
)
FS.time.mat[i,3] <- time.tmp$Elapsed_Time_sec
total.memory.usage[i,"HVGdisp"] <- time.tmp$Total_RAM_Used_MiB
peak.memory.usage[i,"HVGdisp"] <- time.tmp$Peak_RAM_Used_MiB
FS.gene.list[[3]][[i]] <- VariableFeatures(seurat.SD)

# DUBStepR --------------
time.tmp <- peakRAM(
dub.list <- DUBStepR(seurat.SD@assays$RNA@data)
)
FS.time.mat[i,4] <- time.tmp$Elapsed_Time_sec
total.memory.usage[i,"DUBStepR"] <- time.tmp$Total_RAM_Used_MiB
peak.memory.usage[i,"DUBStepR"] <- time.tmp$Peak_RAM_Used_MiB

dub <- dub.list[["optimal.feature.genes"]]
FS.gene.list[[4]][[i]] <- dub
rm(dub,dub.list)

# devianceFS --------------
library(scry)
time.tmp <- peakRAM(
devianceFS_out <- scry::devianceFeatureSelection(object = seurat.SD@assays$RNA@counts)
)
FS.time.mat[i,6] <- time.tmp$Elapsed_Time_sec
total.memory.usage[i,"devianceFS"] <- time.tmp$Total_RAM_Used_MiB
peak.memory.usage[i,"devianceFS"] <- time.tmp$Peak_RAM_Used_MiB

names(devianceFS_out) <- rownames(seurat.SD)
FS.gene.list[[6]][[i]] <- names(sort(devianceFS_out, decreasing = TRUE))
print("devianceFS Finish!")
rm(devianceFS_out)

# TrendVar --------------
## Codes are modified from https://github.com/prabhakarlab/DUBStepR
source("trendVar.R")
library(SingleCellExperiment)
trendVarFS <- function(counts, data) {
  st <- system.time({
    sce <- SingleCellExperiment(list(counts = counts, logcounts = data))
    mgvar <- scran::modelGeneVar(x = sce)
    top.hvgs <- scran::getTopHVGs(mgvar, n = nrow(mgvar))
  })
  
  return(list("var.out" = mgvar, "genes" = top.hvgs, "st" = st))
}
time.tmp <- peakRAM(
trendvar <- trendVarFS(seurat.SD@assays$RNA@counts,seurat.SD@assays$RNA@data)
)

FS.time.mat[i,7] <- time.tmp$Elapsed_Time_sec
total.memory.usage[i,"TrendVar"] <- time.tmp$Total_RAM_Used_MiB
peak.memory.usage[i,"TrendVar"] <- time.tmp$Peak_RAM_Used_MiB
# sum(trendvar[["var.out"]]@listData$FDR < 0.05,na.rm = T)
FS.gene.list[[7]][[i]] <- trendvar["genes"][[1]]
print("TrendVar Finish")
rm(trendvar)

# M3DropDANB --------------
## Codes are modified from https://github.com/prabhakarlab/DUBStepR
source("UseM3D.R")
tryCatch({
  withTimeout({
    time.tmp <- peakRAM(
    M3DropDANB <- UseM3D(seurat.SD@assays$RNA@data,"DANB")
    )
    FS.time.mat[i,8] <- time.tmp$Elapsed_Time_sec
    total.memory.usage[i,"M3DropDANB"] <- time.tmp$Total_RAM_Used_MiB
    peak.memory.usage[i,"M3DropDANB"] <- time.tmp$Peak_RAM_Used_MiB
    
    FS.gene.list[[8]][[i]] <- M3DropDANB$genes
  },timeout = 15000,onTimeout = "error")
}, error = function(e){
  cat("M3DropDANB Error! Simulation ",i,".\n")
  cat("ERROR :",conditionMessage(e),"\n")
})
print("M3DropDANB Finish!")
rm(M3DropDANB)
gc(verbose = F)
