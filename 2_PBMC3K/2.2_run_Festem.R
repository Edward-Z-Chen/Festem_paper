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
set.seed(321)
cl <- makeCluster(getOption("cl.cores", 12))
clusterSetRNGStream(cl, iseed = 321)

pbmc <- readRDS("./results/pbmc3k.rds")
load("./results/pbmc3k_preclustering.RData")
levels(cluster.labels) <- c(1:8,8)
counts <- pbmc@assays$RNA@counts
raw.lib <- pbmc@meta.data$nCount_RNA
counts <- as.matrix(counts)

# For prior with 7 groups
# load("./results/pbmc3k_preclustering_7g.RData")
# For prior with 10 groups
# load("./results/pbmc3k_preclustering_10g.RData")
# Also change the name of the output with corresponding suffix.

# Preprocessing -----------------------------------------------------------

for (i in 1:nrow(counts)){
  tmp <- counts[i,]>=quantile(counts[i,],0.95)
  Q13 <- quantile(counts[i,tmp],c(0.25,0.75))
  if (Q13[2]==Q13[1]) Q13[2] <- Q13[1]+1
  upper <- max(Q13[2]*3-2*Q13[1],1)
  lower <- max(Q13[1]*3-2*Q13[2],0)
  outlier.index <- tmp & (counts[i,]>upper | counts[i,]<lower)
  mean.noout <- mean(counts[i,tmp & (!outlier.index)])
  counts[i,outlier.index] <- round(mean.noout)
}

## Sub-sampling
library.size <- calcNormFactors(counts,lib.size = raw.lib)
sub.sample <- function(x,library.size){
  sample.count <- function(a,library.size){
    rpois(1,a/library.size)
  }
  apply(matrix(x), 1, sample.count,library.size = library.size)
}

for (i in 1:ncol(counts)){
  counts[,i] <- sub.sample(counts[,i],library.size[i])
}

## Outlier2
for (i in 1:nrow(counts)){
  tmp <- counts[i,]>=quantile(counts[i,],0.85)
  if (sum(tmp)!=0){
    Q13 <- quantile(counts[i,tmp],c(0.25,0.75))
    if (Q13[2]==Q13[1]) Q13[2] <- Q13[1]+1
    upper <- max(Q13[2]*3-2*Q13[1],1)
    lower <- max(Q13[1]*3-2*Q13[2],0)
    outlier.index <- tmp & (counts[i,]>upper | counts[i,]<lower)
    mean.noout <- mean(counts[i,tmp & (!outlier.index)])
    counts[i,outlier.index] <- round(mean.noout)
  }
}

# Only consider genes with non-zero expression in at least 30 cells
nonzeros.num <- function(x){sum(x!=0)}
tmp <- apply(counts, 1, nonzeros.num)
counts <- counts[tmp>=30,]


# Run Festem --------------------------------------------------------------
alpha.label <- numeric(nlevels(cluster.labels)-1)
for (g in 1:length(alpha.label)) {
  alpha.label[g] <- sum(cluster.labels==g)/ncol(counts)
}

time.tmp <- Sys.time()
em.result <- parApply(cl,counts,1,em.stat,alpha.ini=rbind(alpha.label,rep(1/8,7)),k0=100,C=1e-3,labels = cluster.labels,group.num = 8,prior.weight=0.05,earlystop = 1e-5)
em.result.9 <- parApply(cl,counts,1,em.stat,alpha.ini=rbind(alpha.label,rep(1/8,7)),k0=100,C=1e-3,labels = cluster.labels,group.num = 8,prior.weight=0.9,earlystop = 1e-5)
print(paste0("Time: ",difftime(Sys.time(),time.tmp,units = "secs")))
save(em.result,em.result.9,file = "./results/pbmc_Festem.RData")

## \gamma = 0.01
em.result <- parApply(cl,counts,1,em.stat,alpha.ini=rbind(alpha.label,rep(1/8,7)),k0=100,C=1e-3,labels = cluster.labels,group.num = 8,prior.weight=0.01,earlystop = 1e-5)
save(em.result,em.result.9,file = "./results/pbmc_Festem_gamma0.01.RData")
stopCluster(cl)