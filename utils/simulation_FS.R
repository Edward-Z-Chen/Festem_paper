library(Seurat)
library(DUBStepR)

all.genes <- rownames(seurat.SD)
## Festem_2g
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

## Festem_3g
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

## Festem_4g
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

## Festem_5g
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

## HVGvst2000
time.tmp <- peakRAM(
seurat.SD <- FindVariableFeatures(seurat.SD, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
)
FS.time.mat[i,2] <- time.tmp$Elapsed_Time_sec
total.memory.usage[i,"HVGvst"] <- time.tmp$Total_RAM_Used_MiB
peak.memory.usage[i,"HVGvst"] <- time.tmp$Peak_RAM_Used_MiB

FS.gene.list[[2]][[i]] <- VariableFeatures(seurat.SD)

## HVGDisp2000
time.tmp <- peakRAM(
seurat.SD <- FindVariableFeatures(object = seurat.SD, selection.method = "disp", nfeatures = 2000,verbose = FALSE)
)
FS.time.mat[i,3] <- time.tmp$Elapsed_Time_sec
total.memory.usage[i,"HVGdisp"] <- time.tmp$Total_RAM_Used_MiB
peak.memory.usage[i,"HVGdisp"] <- time.tmp$Peak_RAM_Used_MiB
FS.gene.list[[3]][[i]] <- VariableFeatures(seurat.SD)

## DUBStepR
time.tmp <- peakRAM(
dub.list <- DUBStepR(seurat.SD@assays$RNA@data)
)
FS.time.mat[i,4] <- time.tmp$Elapsed_Time_sec
total.memory.usage[i,"DUBStepR"] <- time.tmp$Total_RAM_Used_MiB
peak.memory.usage[i,"DUBStepR"] <- time.tmp$Peak_RAM_Used_MiB

dub <- dub.list[["optimal.feature.genes"]]
FS.gene.list[[4]][[i]] <- dub
rm(dub,dub.list)

## HLG
## Adapted from DUBStepR
library(irlba)
irlba_pca_fs <- function (expr_mat, pcs = c(2, 3)) 
{
  ## Taken from DUBStepR
  lognorm <- expr_mat
  nz_genes <- which(Matrix::rowSums(lognorm) != 0)
  lognorm[nz_genes, ] <- lognorm[nz_genes, ]/log(2)
  gene_names <- rownames(lognorm)
  lognorm <- Matrix::Matrix(lognorm, sparse = TRUE)
  rownames(lognorm) <- gene_names
  nc <- ncol(lognorm)
  expression_means <- Matrix::rowMeans(lognorm)
  expression_vars <- Matrix::rowMeans((lognorm - expression_means)^2) * 
    (nc/(nc - 1))
  genes_to_keep <- expression_vars > 0
  lognorm <- lognorm[genes_to_keep, ]
  expression_means <- expression_means[genes_to_keep]
  expression_vars <- expression_vars[genes_to_keep]
  irlba_pca_res <- irlba::irlba(Matrix::t(lognorm), nu = 0, center = expression_means, nv = 30,
                                scale = sqrt(expression_vars), right_only = TRUE)$v
  row.names(irlba_pca_res) <- row.names(lognorm)
  if (length(pcs) > 1) {
    score <- Matrix::rowSums(abs(irlba_pca_res[, pcs]))
  }
  else {
    score <- abs(irlba_pca_res[, pcs])
  }
  names(score) = gene_names[genes_to_keep]
  return(sort(-score))
}

time.tmp <- peakRAM(
PCA_Score <- irlba_pca_fs(expr_mat = seurat.SD@assays$RNA@data, pcs = 1:5)
)
FS.time.mat[i,5] <- time.tmp$Elapsed_Time_sec
total.memory.usage[i,"HLG_5pc"] <- time.tmp$Total_RAM_Used_MiB
peak.memory.usage[i,"HLG_5pc"] <- time.tmp$Peak_RAM_Used_MiB

FS.gene.list[[5]][[i]] <- names(PCA_Score)

time.tmp <- peakRAM(
  PCA_Score <- irlba_pca_fs(expr_mat = seurat.SD@assays$RNA@data, pcs = 1:10)
)
FS.time.mat[i,9] <- time.tmp$Elapsed_Time_sec
total.memory.usage[i,"HLG_10pc"] <- time.tmp$Total_RAM_Used_MiB
peak.memory.usage[i,"HLG_10pc"] <- time.tmp$Peak_RAM_Used_MiB

FS.gene.list[[9]][[i]] <- names(PCA_Score)

print("HLG Finish!")
rm(PCA_Score)

##devianceFS
## This can enable batch!!! (in Shiyang)
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

## TrendVar
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

##M3DropDANB
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