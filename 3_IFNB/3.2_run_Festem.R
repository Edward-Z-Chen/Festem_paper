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

ifnb <- readRDS("./results/ifnb_ctrl.rds")
### This prior is generated with all genes, 20pca, 0.7reso on my laptop (code in 3.1)
load("./results/ifnb_ctrl_preclustering.RData")

### This prior is generated with all genes, 20pca, 0.6reso on my laptop (code in 3.1)
# load("./results/ifnb_ctrl_preclustering_13g.RData")

### This prior is generated with all genes, 20pca, 1reso on my laptop (code in 3.1)
# load("./results/ifnb_ctrl_preclustering_17g.RData")
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

save(em.result,em.result.9,file = "./results/ifnb_ctrl_Festem.RData")

## \gamma = 0.01
time.tmp <- Sys.time()
em.result <- parApply(cl,counts,1,em.stat,alpha.ini=rbind(alpha.label,rep(1/nlevels(cluster.labels),nlevels(cluster.labels)-1)),k0=100,C=1e-3,labels = cluster.labels,group.num = nlevels(cluster.labels),prior.weight=0.01,earlystop = 1e-5)
print(paste0("Time cost: ",difftime(Sys.time(),time.tmp,units = "secs")))
save(em.result,em.result.9,file = "./results/ifnb_ctrl_Festem_gamma0.01.RData")

stopCluster(cl)