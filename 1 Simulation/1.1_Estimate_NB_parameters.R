# Estimate parameters of NB distribution from PBMC3K dataset --------------
library(BiocParallel)
library(BiocGenerics)
library(parallel)
library(MASS)
library(nloptr)
library(Seurat)

pbmc <- readRDS("pbmc_tutorial.rds")
counts <- pbmc@assays$RNA@counts
cluster.label <- pbmc@active.ident
rm(pbmc)

counts <- as.matrix(counts)
counts0 <- counts[,cluster.label==levels(cluster.label)[1]]
counts1 <- counts[,cluster.label==levels(cluster.label)[2]]
counts2 <- counts[,cluster.label==levels(cluster.label)[3]]
counts3 <- counts[,cluster.label==levels(cluster.label)[4]]
counts4 <- counts[,cluster.label==levels(cluster.label)[5]]
counts5 <- counts[,cluster.label==levels(cluster.label)[6]]
counts6 <- counts[,cluster.label==levels(cluster.label)[7]]
rm(counts)

fit.nb <- function(x){
  require(MASS)
  require(nloptr)
  # this function returns the MLE of mean and r
  obj.f <- function(theta){-sum(dnbinom(x,theta[2],mu = theta[1],log = T))}
  deriv.pl <- function(theta){nl.grad(x0 = theta,fn = obj.f)}
  nloptr(x0 = c(max(mean(x),1e-6),5),eval_f = obj.f,lb=rep(0,2),ub = rep(300,2),
         opts = list("algorithm"="NLOPT_LN_NELDERMEAD","ftol_rel" = 1e-9,
                     "maxeval" = 5000,"maxtime" = 200,"xtol_rel" = 1e-4))$solution
}

cl <- makeCluster(getOption("cl.cores", 12))
counts.nb0 <- parApply(cl,counts0, 1, fit.nb)
counts.nb1 <- parApply(cl,counts1, 1, fit.nb)
counts.nb2 <- parApply(cl,counts2, 1, fit.nb)
counts.nb3 <- parApply(cl,counts3, 1, fit.nb)
counts.nb4 <- parApply(cl,counts4, 1, fit.nb)
counts.nb5 <- parApply(cl,counts5, 1, fit.nb)
counts.nb6 <- parApply(cl,counts6, 1, fit.nb)
stopCluster(cl)
save(counts.nb0,counts.nb1,counts.nb2,counts.nb3,counts.nb4,counts.nb5,counts.nb6,file = "counts_mean_r.RData")
