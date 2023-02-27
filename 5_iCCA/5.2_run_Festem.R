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
library(EMDE)
set.seed(321)
mc.reset.stream()
load("./results/iCCA_preclustering.RData")
nonepi <- readRDS("NonEpi.rds")

counts <- nonepi@assays$RNA@counts
batch.id <- factor(nonepi@meta.data$Patient)
rm(nonepi)
cluster.labels <- factor(cluster.labels)
levels(cluster.labels) <- 1:nlevels(cluster.labels)

all.genes <- rownames(counts)

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
sub.sample <- function(x){
  library.size <- x[1]
  x <- x[-1]
  sample.count <- function(a,library.size){
    rpois(1,a/library.size)
  }
  apply(matrix(x), 1, sample.count,library.size = library.size)
}

cl <- makeCluster(getOption("cl.cores", 12))
clusterSetRNGStream(cl, iseed = 321)
em.result <- vector("list",nlevels(batch.id))
em.result.9 <- vector("list",nlevels(batch.id))
for (B in 1:nlevels(batch.id)){
  print(paste0("Batch: ",levels(batch.id)[B]))
  counts.tmp <- counts[,batch.id==(levels(batch.id)[B])]
  cluster.labels.tmp <- cluster.labels[batch.id==(levels(batch.id)[B])]
  # Removing those clusters with no more than 10 cells in this batch
  cluster.labels.tmp <- factor(cluster.labels.tmp)
  counts.tmp <- counts.tmp[,summary(cluster.labels.tmp)[cluster.labels.tmp]>10]
  cluster.labels.tmp <- cluster.labels.tmp[summary(cluster.labels.tmp)[cluster.labels.tmp]>10]
  cluster.labels.tmp <- factor(cluster.labels.tmp)
  levels(cluster.labels.tmp) <- 1:nlevels(cluster.labels.tmp)
  
  ## Outlier
  counts.tmp <- t(parApply(cl,counts.tmp,1,rm.outlier,percent = 0.95))
  rownames(counts.tmp) <- rownames(counts)
  
  ## Sub-sampling
  library.size <- calcNormFactors(counts.tmp)
  counts.tmp <- parApply(cl,rbind(library.size,counts.tmp),2,sub.sample)
  counts.tmp <- t(parApply(cl,counts.tmp,1,rm.outlier,percent = 0.90))
  rownames(counts.tmp) <- rownames(counts)
  
  nonzeros.num <- function(x){sum(x!=0)}
  tmp <- apply(counts.tmp, 1, nonzeros.num)
  min.cell <- min(0.01*ncol(counts.tmp),30)
  sum(tmp>=min.cell)
  counts.tmp <- counts.tmp[tmp>=min.cell,]
  
  alpha.label <- numeric(nlevels(cluster.labels.tmp)-1)
  for (g in 1:length(alpha.label)) {
    alpha.label[g] <- sum(cluster.labels.tmp==g)/ncol(counts.tmp)
  }
  
  time.tmp <- Sys.time()
  em.result[[B]] <- parApply(cl,counts.tmp,1,em.stat,alpha.ini=rbind(alpha.label,rep(1/nlevels(cluster.labels.tmp),length(alpha.label))),k0=100,C=1e-3,labels = cluster.labels.tmp,
                             group.num = nlevels(cluster.labels.tmp),prior.weight=0.05,earlystop = 1e-4)
  print(paste0("Time cost: ",difftime(Sys.time(),time.tmp,units = "secs")))
  
  time.tmp <- Sys.time()
  em.result.9[[B]] <- parApply(cl,counts.tmp,1,em.stat,alpha.ini=rbind(alpha.label,rep(1/nlevels(cluster.labels.tmp),length(alpha.label))),k0=100,C=1e-3,labels = cluster.labels.tmp,
                               group.num = nlevels(cluster.labels.tmp),prior.weight=0.9,earlystop = 1e-4)
  print(paste0("Time cost: ",difftime(Sys.time(),time.tmp,units = "secs")))
  
  save(all.genes,em.result,em.result.9,file = "./results/iCCA_Festem.RData")
}
save(all.genes,em.result,em.result.9,file = "./results/iCCA_Festem.RData")
