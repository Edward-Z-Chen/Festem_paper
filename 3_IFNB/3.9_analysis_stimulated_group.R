# Download preprocessed data -----------------------------------------------------------
library(Seurat)
library(SeuratData)
InstallData("ifnb")
data("ifnb")
ifnb <- subset(ifnb,subset = orig.ident!="IMMUNE_CTRL")
saveRDS(ifnb,file = "./results/ifnb_stim.rds")


# Preclustering -----------------------------------------------------------
ifnb <- readRDS("./results/ifnb_stim.rds")
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb,nfeatures = 8000)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb,verbose = F)
ifnb <- FindNeighbors(ifnb, dims = 1:20)
ifnb <- FindClusters(ifnb, resolution = 0.7)
ifnb <- RunTSNE(ifnb,dims = 1:20)
DimPlot(ifnb, reduction = "tsne",label = T)

cluster.labels <- ifnb@meta.data$seurat_clusters
levels(cluster.labels) <- 1:nlevels(cluster.labels)


# Festem ------------------------------------------------------------------

library(SingleCellExperiment)
library(scater)
library(stringr)
library(stats)
library(BiocParallel)
library(BiocGenerics)
library(parallel)
library(MASS)
library(nloptr)
library(edgeR)
library(EMDE)
library(Seurat)

set.seed(321)
cl <- makeCluster(getOption("cl.cores", 12))
clusterSetRNGStream(cl, iseed = 321)

ifnb <- readRDS("./results/ifnb_stim.rds")
levels(cluster.labels) <- 1:nlevels(cluster.labels)
counts <- ifnb@assays$RNA@counts
raw.lib <- ifnb@meta.data$nCount_RNA
counts <- as.matrix(counts)

## Outlier
rm.outlier <- function(x,percent){
  tmp <- (x>=quantile(x,percent))
  Q13 <- quantile(x[tmp],c(0.25,0.75))
  if (Q13[2]==Q13[1]) Q13[2] <- Q13[1]+1
  upper <- max(Q13[2]*3-2*Q13[1],1)
  lower <- max(Q13[1]*3-2*Q13[2],0)
  outlier.index <- tmp & (x>upper | x<lower)
  mean.noout <- mean(x[tmp & (!outlier.index)])
  x[outlier.index] <- round(mean.noout)
  cat(sum(outlier.index),"\n")
  x
}
counts <- t(parApply(cl,counts,1,rm.outlier,percent = 0.95))

## Sub-sampling
library.size <- calcNormFactors(counts,lib.size = raw.lib)
sub.sample <- function(x){
  library.size <- x[1]
  x <- x[-1]
  sample.count <- function(a,library.size){
    rpois(1,a/library.size)
  }
  apply(matrix(x), 1, sample.count,library.size = library.size)
}

counts <- parApply(cl,rbind(library.size,counts),2,sub.sample)
counts <- t(parApply(cl,counts,1,rm.outlier,percent = 0.85))
rownames(counts) <- rownames(ifnb)

nonzeros.num <- function(x){sum(x!=0)}
tmp <- apply(counts, 1, nonzeros.num)
sum(tmp>=30)
counts <- counts[tmp>=30,]

alpha.label <- numeric(nlevels(cluster.labels)-1)
for (i in 1:(nlevels(cluster.labels)-1)) {
  alpha.label[i] <- sum(cluster.labels==i)/ncol(counts)
}


time.tmp <- Sys.time()
em.result <- parApply(cl,counts,1,em.stat,alpha.ini=rbind(alpha.label,rep(1/nlevels(cluster.labels),nlevels(cluster.labels)-1)),k0=100,C=1e-3,labels = cluster.labels,group.num = nlevels(cluster.labels),prior.weight=0.05,earlystop = 1e-5)
print(paste0("Time cost: ",difftime(Sys.time(),time.tmp,units = "secs")))

time.tmp <- Sys.time()
em.result.9 <- parApply(cl,counts,1,em.stat,alpha.ini=rbind(alpha.label,rep(1/nlevels(cluster.labels),nlevels(cluster.labels)-1)),k0=100,C=1e-3,labels = cluster.labels,group.num = nlevels(cluster.labels),prior.weight=0.9,earlystop = 1e-5)
print(paste0("Time cost: ",difftime(Sys.time(),time.tmp,units = "secs")))

save(em.result,em.result.9,file = "./results/ifnb_stim_Festem.RData")

# Clustering and TSNE ------------------------------------------------------------------
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

ifnb <- readRDS("./results/ifnb_stim.rds")
ifnb <- NormalizeData(ifnb)
ifnb <- ifnb[,!ifnb@meta.data$seurat_annotations%in%c("Eryth","Mk")]

ifnb <- ScaleData(ifnb,features = EM[1:2500])
ifnb <- RunPCA(ifnb, verbose = FALSE,features = EM[1:2500])
ifnb <- FindNeighbors(object = ifnb, dims = 1:25)
ifnb <- FindClusters(object = ifnb, resolution = 1.5)
ifnb <- RunTSNE(ifnb, reduction = "pca", dims = 1:25)
label <- ifnb@active.ident
tsne <- ifnb@reductions[["tsne"]]
save(label,tsne,file = "./results/ifnb_stim_clustering_tSNE.RData")
