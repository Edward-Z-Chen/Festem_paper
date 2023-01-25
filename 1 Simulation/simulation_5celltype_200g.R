library(SingleCellExperiment)
library(scater)
library(stringr)
library(stats)
library(BiocParallel)
library(BiocGenerics)
library(parallel)
library(MASS)
library(nloptr)
library(Seurat)
library(R.utils)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(edgeR)
library(EMDE)
library(peakRAM)
library(pracma)
library(lcmix)
library(Matrix)

load("counts_mean_r.RData")

ecdf.em <- function(x){
  pchisq(x,2)
}

my.wilcox <- function(x,label){
  # labels should be T or F
  wilcox.test(x[label],x[!label])$p.value
}

umap.for.plot <- function(umap,cluster){
  umap <- umap@cell.embeddings
  umap <- as.data.frame(umap)
  cbind(umap,cluster = cluster)
}

ratio <- 0.5
DE.gene <- 200
nonDE.gene <- 19800
samples <- 3000
cor_num <- 200
rho <- 0.99
corr <- matrix(rho, nrow = cor_num, ncol = cor_num)
corr <- corr + diag(rep(1-rho, cor_num))
mu.list <- numeric(DE.gene+nonDE.gene)
sigma.list <- numeric(DE.gene+nonDE.gene)
sigma.list2 <- numeric(DE.gene)
counts.mean <- counts.nb0[1,]
counts.r <- counts.nb0[2,]
counts.mean2 <- counts.nb1[1,]
counts.r2 <- counts.nb1[2,]
seedlist <- c(234,654,953,1245,563,11143,632,920,804,336,
              154,4373,6213,791,8832,9012,1255,7864,15627,2341)
B <- 20
save.image("./Simulation_V3/NB_200DE_5type_workspace.RData")

### Adjusted p-value list
adjpvalue.list <- vector("list",21)
names(adjpvalue.list) <- c("Festem_0.05","Festem_0.05_f",
                           "DESeq2","DESeq2_full",
                           "EdgeR","EdgeR_full",
                           "MAST","Wilcoxon",
                           "ZIAQ","scDD","sigEMD","SCDE",
                           "ROSeq","Festem_3group_0.05","Festem_3group_0.05f",
                           "Festem_7group_0.05","Festem_7group_0.05f",
                           "Festem_9group_0.05","Festem_9group_0.05f",
                           "singleCellHaystack_PCA",
                           "singleCellHaystack_UMAP")
for (i in 1:length(adjpvalue.list)){
  adjpvalue.list[[i]] <- matrix(nrow = B, ncol = DE.gene+nonDE.gene)
}
time.mat <- matrix(nrow = B, ncol = length(adjpvalue.list))
colnames(time.mat) <- names(adjpvalue.list)

cluster.acc <- numeric(B)
cluster.labels.mat <- matrix(nrow = B, ncol = samples)

FS.gene.list <- vector("list",12)
names(FS.gene.list) <- c("Festem_5","HVGvst","HVGdisp","DUBStepR",
                         "HLG","devianceFS","TrendVar","M3DropDANB",
                         "DUBStepR_full","Festem_3","Festem_7","Festem_9")
for (i in 1:length(FS.gene.list)){
  FS.gene.list[[i]] <- vector("list",B)
}
FS.time.mat <- matrix(nrow = B,ncol = length(FS.gene.list))
colnames(FS.time.mat) <- names(FS.gene.list)


FS.label.list <- vector("list",13)
names(FS.label.list) <- c("Festem","HVGvst","HVGdisp","DUBStepR",
                          "HLG","devianceFS","TrendVar","M3DropDANB",
                          "DUBStepR_full","NO FS","Festem_3","Festem_7","Festem_9")
for (i in 1:length(FS.label.list)){
  FS.label.list[[i]] <- matrix(nrow = B, ncol = samples)
}

FS.ARI.frame <- matrix(nrow = B,ncol = length(FS.label.list))
colnames(FS.ARI.frame) <- names(FS.label.list)

FS.label.kmeans.list <- vector("list",13)
names(FS.label.kmeans.list) <- c("Festem","HVGvst","HVGdisp","DUBStepR",
                                 "HLG_5pc","devianceFS","TrendVar","M3DropDANB",
                                 "HLG_10pc","NO FS","Festem_3","Festem_4","Festem_5")
for (i in 1:length(FS.label.kmeans.list)){
  FS.label.kmeans.list[[i]] <- matrix(nrow = B, ncol = samples)
}

FS.ARI.kmeans.frame <- matrix(nrow = B,ncol = length(FS.label.kmeans.list))
colnames(FS.ARI.kmeans.frame) <- names(FS.label.kmeans.list)

FS.SI.frame <- matrix(nrow = B,ncol = length(FS.label.list))
colnames(FS.SI.frame) <- names(FS.label.list)

FS.DB.frame <- matrix(nrow = B,ncol = length(FS.label.list))
colnames(FS.DB.frame) <- names(FS.label.list)

FS.label.true.list <- vector("list",13)
names(FS.label.true.list) <- c("Festem","HVGvst","HVGdisp","DUBStepR",
                               "HLG_5pc","devianceFS","TrendVar","M3DropDANB",
                               "HLG_10pc","NO FS","Festem_3","Festem_4","Festem_5")
for (i in 1:length(FS.label.true.list)){
  FS.label.true.list[[i]] <- matrix(nrow = B, ncol = samples)
}

FS.ARI.true.frame <- matrix(nrow = B,ncol = length(FS.label.true.list))
colnames(FS.ARI.true.frame) <- names(FS.label.true.list)

FS.label.kmeans.true.list <- vector("list",13)
names(FS.label.kmeans.true.list) <- c("Festem","HVGvst","HVGdisp","DUBStepR",
                                      "HLG_5pc","devianceFS","TrendVar","M3DropDANB",
                                      "HLG_10pc","NO FS","Festem_3","Festem_4","Festem_5")
for (i in 1:length(FS.label.kmeans.true.list)){
  FS.label.kmeans.true.list[[i]] <- matrix(nrow = B, ncol = samples)
}

FS.ARI.kmeans.true.frame <- matrix(nrow = B,ncol = length(FS.label.kmeans.true.list))
colnames(FS.ARI.kmeans.true.frame) <- names(FS.label.kmeans.true.list)

FS.SI.true.frame <- matrix(nrow = B,ncol = length(FS.label.true.list))
colnames(FS.SI.true.frame) <- names(FS.label.true.list)

FS.DB.true.frame <- matrix(nrow = B,ncol = length(FS.label.true.list))
colnames(FS.DB.true.frame) <- names(FS.label.true.list)

all.method <- c("Festem_0.05","Festem_0.05f",
                "DESeq2","DESeq2_full",
                "EdgeR","EdgeR_full",
                "MAST","Wilcoxon",
                "ZIAQ","scDD","sigEMD","SCDE",
                "ROSeq","Festem_3group_0.05","Festem_3group_0.05f",
                "Festem_7group_0.05","Festem_7group_0.05f",
                "Festem_9group_0.05","Festem_9group_0.05f",
                "singleCellHaystack_PCA",
                "singleCellHaystack_UMAP",
                "HVGvst","HVGdisp","DUBStepR",
                "HLG_5pc","devianceFS","TrendVar","M3DropDANB","HLG_10pc")
total.memory.usage <- matrix(nrow = B,ncol = length(all.method))
peak.memory.usage <- matrix(nrow = B,ncol = length(all.method))
colnames(total.memory.usage) <- all.method
colnames(peak.memory.usage) <- all.method

for (i in 1:B){
  ### Unique dispersion
  print(paste0("Round ",i,":"))
  set.seed(seedlist[i])
  cl <- makeCluster(getOption("cl.cores", 12))
  clusterSetRNGStream(cl, iseed = seedlist[i])
  counts <- matrix(rep(0,(DE.gene+nonDE.gene)*samples),ncol = samples)
  colnames(counts) <- paste("Cell",seq(1:samples),sep = "")
  rownames(counts) <- paste("Gene",seq(1:(DE.gene+nonDE.gene)),sep = "")
  ## non-DE genes
  generate.counts <- function(mu,n,sig){
    rnbinom(n,sig,mu=mu)
  }
  ### Correlated Part
  mu <- sample(counts.mean[counts.mean>= 0.5 & counts.mean<=30],cor_num,replace = T)
  dispersion <- sample(counts.r[counts.r>=0.2 & counts.r<=50],cor_num,replace = T)
  gamma_dta <- rmvgamma(samples, shape = dispersion, rate = dispersion/(mu), corr = corr)
  counts[1:cor_num,] <- (apply(gamma_dta, 1, function(lam){rpois(length(lam), lambda = lam)}))
  mu.list[1:cor_num] <- mu
  sigma.list[1:cor_num] <- dispersion
  
  mu <- sample(counts.mean[counts.mean>=0.5 & counts.mean<=30],nonDE.gene - cor_num,replace = T)
  mu.list[(cor_num+1):nonDE.gene] <- mu
  sigma.list[(cor_num+1):nonDE.gene] <- sample(counts.r[counts.r>=0.2 & counts.r<=50],nonDE.gene - cor_num,replace = T)
  for (k in (cor_num+1):nonDE.gene){
    counts[k,] <- rnbinom(samples,sigma.list[k],mu = mu.list[k])
  }
  
  ## DE genes
  mu <- sample(counts.mean[counts.mean>=0.5 & counts.mean<=30],DE.gene,replace = T)
  log.foldchange <- runif(DE.gene,1.5,2.5)
  
  ## Only considers up-regulation
  mu2 <- mu*2^log.foldchange
  mu.list[(nonDE.gene+1):(nonDE.gene+DE.gene)] <- mu
  sigma.list[(nonDE.gene+1):(nonDE.gene+DE.gene)] <- sample(counts.r[counts.r>=0.2 & counts.r<=20],DE.gene,replace = T)
  sigma.list2 <- sample(counts.r2[counts.r2>=0.2 & counts.r<=20],DE.gene,replace = T)
  for (k in 1:DE.gene){
    if (k!=1){
      counts[(k+nonDE.gene),1:(samples*ratio)] <- rnbinom((samples*ratio),sigma.list[k+nonDE.gene],mu = mu[k])
      if ((k %% 4) == 0) index.tmp <- 1501:1875
      if ((k %% 4) == 1) index.tmp <- 1876:2250
      if ((k %% 4) == 2) index.tmp <- 2251:2625
      if ((k %% 4) == 3) index.tmp <- 2626:3000
      counts[(k+nonDE.gene),index.tmp] <- rnbinom((samples*(1-ratio)*0.25),sigma.list[k+nonDE.gene],mu = mu[k])
      counts[(k+nonDE.gene),((samples*ratio+1):samples)[!((samples*ratio+1):samples %in% index.tmp)]] <- rnbinom((samples*(1-ratio)*0.75),sigma.list2[k],mu = mu2[k])
    }
    else {
      counts[(k+nonDE.gene),1:(samples*ratio)] <- rnbinom((samples*ratio),sigma.list[k+nonDE.gene],mu = mu[k])
      counts[(k+nonDE.gene),(samples*ratio+1):samples] <- rnbinom((samples*(1-ratio)),sigma.list2[k],mu = mu2[k])
    }
  }
  save(counts,mu.list,log.foldchange,sigma.list,sigma.list2,file = paste0("./Simulation_V3/NB_200DE_5type_",i,".RData"))
  print("Finish generating counts!")
  
  # Clustering
  seurat.SD <- CreateSeuratObject(counts)
  seurat.SD <- NormalizeData(seurat.SD)
  
  SD <- SingleCellExperiment(
    assays = list(counts = counts)
  )
  SD <- logNormCounts(SD)
  assayNames(SD)[2] <- "normcounts"
  
  seurat.SD <- FindVariableFeatures(seurat.SD,selection.method = "vst",nfeatures = 1000)
  seurat.SD <- ScaleData(seurat.SD)
  seurat.SD <- RunPCA(seurat.SD, verbose = FALSE)
  seurat.SD <- FindNeighbors(object = seurat.SD, dims = 1:10,reduction = "pca")
  
  for (k in 1:100){
    seurat.SD <- FindClusters(object = seurat.SD, resolution = 0.01*k, verbose = FALSE)
    if (nlevels(seurat.SD@active.ident)>=5){
      print(paste0("Resolution: ", 0.01*k,"; Group number: ",nlevels(seurat.SD@active.ident)))
      break
    }
  }
  
  # hvg <- VariableFeatures(seurat.SD)
  # clustering <- kmeans(t(seurat.SD@assays$RNA@data[hvg,]),5)
  seurat.SD <- AddMetaData(seurat.SD,seurat.SD@active.ident,"Louvain")
  seurat.SD <- AddMetaData(seurat.SD,c(rep(1,samples * ratio),rep("2_1",375),rep("2_2",375),
                                       rep("2_3",375),rep("2_4",375)),"truth")
  
  colData(SD)[,1] <- as.numeric(seurat.SD@active.ident)
  names(SD@colData@listData)[1] <- "condition"
  
  cluster.labels.mat[i,] <- as.numeric(seurat.SD@active.ident)
  cluster.labels <- factor(as.numeric(seurat.SD@active.ident))
  cluster.acc[i] <- mclust::adjustedRandIndex(as.numeric(cluster.labels),
                                              seurat.SD@meta.data$truth)
  print("Finish clustering!")
  
  ############### UMAP ##################
  seurat.SD <- FindVariableFeatures(seurat.SD,selection.method = "vst",nfeatures = 1000)
  seurat.SD <- ScaleData(seurat.SD)
  seurat.SD <- RunPCA(seurat.SD)
  seurat.SD <- RunUMAP(seurat.SD,dims = 1:10)
  if (i == B | i == 1){
    pdf("./Simulation_V3/NB_200DE_5type_UMAP_initial.pdf",width=4,height=4,onefile = T)
    print(DimPlot(seurat.SD,label = T, group.by = "truth",pt.size = 1))
    print(DimPlot(seurat.SD,label = T, group.by = "Louvain",pt.size = 1))
    dev.off()
  }
  
  ################## SingleCellHaystack ########################
  library(singleCellHaystack)
  tryCatch({
    withTimeout({
      time.tmp <- peakRAM(
        haystack_umap <- haystack(seurat.SD,coord = "umap",method = "2D")
      )
      time.mat[i,"singleCellHaystack_UMAP"] <- time.tmp$Elapsed_Time_sec
      total.memory.usage[i,"singleCellHaystack_UMAP"] <- time.tmp$Total_RAM_Used_MiB
      peak.memory.usage[i,"singleCellHaystack_UMAP"] <- time.tmp$Peak_RAM_Used_MiB
      adjpvalue.list[["singleCellHaystack_UMAP"]][i,] <- exp(haystack_umap$results$log.p.adj)
      
      time.tmp <- peakRAM(
        haystack_pca <- haystack(seurat.SD,coord = "pca",dims = 1:10,method = "highD")
      )
      time.mat[i,"singleCellHaystack_PCA"] <- time.tmp$Elapsed_Time_sec
      total.memory.usage[i,"singleCellHaystack_PCA"] <- time.tmp$Total_RAM_Used_MiB
      peak.memory.usage[i,"singleCellHaystack_PCA"] <- time.tmp$Peak_RAM_Used_MiB
      adjpvalue.list[["singleCellHaystack_PCA"]][i,] <- exp(haystack_pca$results$log.p.adj)
      print("singleCellHaystack Finish!")
    },timeout = 30000,onTimeout = "error")
  }, error = function(e){
    cat("singleCellHaystack Error! Simulation ",i,".\n")
    cat("ERROR :",conditionMessage(e),"\n")
  })
  save(adjpvalue.list,time.mat,cluster.acc,cluster.labels.mat,total.memory.usage,peak.memory.usage,file = "./Simulation_V3/NB_200DE_5type_DEG.RData")
  rm(haystack_umap,haystack_pca)
  gc(verbose = F)
  
  
  ################## ROSeq ########################
  library(ROSeq)
  library(limma)
  library(edgeR)
  tryCatch({
    withTimeout({
      time.tmp <- peakRAM({
        counts2<-limma::voom(ROSeq::TMMnormalization(counts))
        
        roseq.tmp <- matrix(nrow = nlevels(cluster.labels),ncol = nrow(counts))
        for (l in 1:nlevels(cluster.labels)){
          condition.tmp <- cluster.labels == l
          condition.tmp <- as.numeric(condition.tmp)+1
          roseq.tmp[l,] <- ROSeq(countData=counts2$E, condition = condition.tmp, numCores=6)[,"pVals"]
        }
      })
      
      time.mat[i,"ROSeq"] <- time.tmp$Elapsed_Time_sec*6
      total.memory.usage[i,"ROSeq"] <- time.tmp$Total_RAM_Used_MiB
      peak.memory.usage[i,"ROSeq"] <- time.tmp$Peak_RAM_Used_MiB
      
      roseq.tmp <- apply(roseq.tmp,2,min,na.rm = T)
      roseq.tmp <- roseq.tmp * nlevels(cluster.labels)
      adjpvalue.list[["ROSeq"]][i,] <- p.adjust(roseq.tmp,"BH")
      print("ROSeq Finish!")
    },timeout = 15000,onTimeout = "error")
  }, error = function(e){
    cat("ROSeq Error! Simulation ",i,".\n")
    cat("ERROR :",conditionMessage(e),"\n")
  })
  save(adjpvalue.list,time.mat,cluster.acc,cluster.labels.mat,total.memory.usage,peak.memory.usage,file = "./Simulation_V3/NB_200DE_5type_DEG.RData")
  rm(counts2,roseq.tmp,condition.tmp)
  gc(verbose = F)
  # detach("package:ROSeq",unload = TRUE)
  
  ## EM test
  alpha.label <- numeric(nlevels(cluster.labels)-1)
  for (g in 1:length(alpha.label)) {
    alpha.label[g] <- sum(cluster.labels==g)/ncol(seurat.SD)
  }
  
  # #### Preprocessing
  # time.tmp <- Sys.time()
  # counts.tmp <- t(parApply(cl,counts,1,rm.outlier,percent = 0.95))
  # rownames(counts.tmp) <- rownames(counts)
  # library.size <- calcNormFactors(counts.tmp)
  # counts.tmp <- parApply(cl,rbind(library.size,counts.tmp),2,sub.sample)
  # # save(counts.tmp,file = paste0("colo_subsample_new",B,".RData"))
  # counts.tmp <- t(parApply(cl,counts.tmp,1,rm.outlier,percent = 0.90))
  # rownames(counts.tmp) <- rownames(counts)
  # pre.time <- difftime(Sys.time(),time.tmp,units = "secs")*12
  
  # ################## EM_0 ########################
  # time.tmp <- Sys.time()
  # em.result <- parApply(cl,seurat.SD@assays$RNA@counts,1,em.stat,
  #                       alpha.ini=rbind(alpha.label,c(rep(1/nlevels(cluster.labels),length(alpha.label)))),
  #                       k0=150,C=1e-3,labels = cluster.labels,group.num = nlevels(cluster.labels),
  #                       prior.weight = 0,earlystop = 1e-5)
  # time.mat[i,1] <- difftime(Sys.time(),time.tmp,units = "secs")*12
  # adjpvalue.list[[1]][i,] <- p.adjust(em.result[1,],"BH")
  
  ################## EM_0.05 ########################
  time.tmp <- peakRAM(
    em.result <- parApply(cl,counts,1,em.stat,
                          alpha.ini=rbind(alpha.label,c(rep(1/nlevels(cluster.labels),length(alpha.label)))),
                          k0=150,C=1e-3,labels = cluster.labels,group.num = nlevels(cluster.labels),
                          prior.weight = 0.05,earlystop = 1e-5)
  )
  time.mat[i,"Festem_0.05"] <- time.tmp$Elapsed_Time_sec*12
  total.memory.usage[i,"Festem_0.05"] <- time.tmp$Total_RAM_Used_MiB
  peak.memory.usage[i,"Festem_0.05"] <- time.tmp$Peak_RAM_Used_MiB
  adjpvalue.list[["Festem_0.05"]][i,] <- p.adjust(em.result[1,],"BH")
  
  ################## EM_0.05f ########################
  time.tmp <- peakRAM(
    em.result.9 <- parApply(cl,counts,1,em.stat,
                            alpha.ini=rbind(alpha.label,c(rep(1/nlevels(cluster.labels),length(alpha.label)))),
                            k0=150,C=1e-3,labels = cluster.labels,group.num = nlevels(cluster.labels),
                            prior.weight = 0.9,earlystop = 1e-5)
  )
  time.mat[i,"Festem_0.05_f"] <- time.tmp$Elapsed_Time_sec*12+time.mat[i,"Festem_0.05"]
  total.memory.usage[i,"Festem_0.05f"] <- max(time.tmp$Total_RAM_Used_MiB,total.memory.usage[i,"Festem_0.05"])
  peak.memory.usage[i,"Festem_0.05f"] <- max(time.tmp$Peak_RAM_Used_MiB,peak.memory.usage[i,"Festem_0.05"])
  
  adjpvalue.list[["Festem_0.05_f"]][i,] <- adjpvalue.list[["Festem_0.05"]][i,]
  adjpvalue.list[["Festem_0.05_f"]][i,em.result.9[2,]<0 & (!is.na(em.result.9[2,])) ] <- NA
  print("Festem Test Finish!")
  save(adjpvalue.list,time.mat,cluster.acc,cluster.labels.mat,total.memory.usage,peak.memory.usage,file = "./Simulation_V3/NB_200DE_5type_DEG.RData")
  
  
  ################## Wilcoxon ########################
  tryCatch({
    withTimeout({
      time.tmp <- peakRAM({
        wil.tmp <- matrix(nrow = nlevels(cluster.labels),ncol = nrow(counts))
        for (l in 1:nlevels(cluster.labels)){
          wil.tmp[l,] <- parApply(cl,counts,1,my.wilcox,label = (cluster.labels==l))
        }
      })
      wil.tmp <- apply(wil.tmp,2,min,na.rm = T)
      wil.tmp <- wil.tmp * nlevels(cluster.labels)
      time.mat[i,"Wilcoxon"] <- time.tmp$Elapsed_Time_sec*12
      total.memory.usage[i,"Wilcoxon"] <- time.tmp$Total_RAM_Used_MiB
      peak.memory.usage[i,"Wilcoxon"] <- time.tmp$Peak_RAM_Used_MiB
      adjpvalue.list[["Wilcoxon"]][i,] <- p.adjust(wil.tmp,"BH")
      print("Wilcoxon Finish!")
    },timeout = 15000,onTimeout = "error")
  }, error = function(e){
    cat("Wilcoxon Error! Simulation ",i,".\n")
    cat("ERROR :",conditionMessage(e),"\n")
  })
  save(adjpvalue.list,time.mat,cluster.acc,cluster.labels.mat,total.memory.usage,peak.memory.usage,file = "./Simulation_V3/NB_200DE_5type_DEG.RData")
  stopCluster(cl)
  rm(wil.tmp)
  gc(verbose = F)
  
  snowparam <- SnowParam(workers = 12, type = "SOCK", log = T)
  register(snowparam, default = TRUE)
  registered()
  
  ################## DESeq2 (filter & full) ########################
  library(DESeq2)
  tryCatch({
    withTimeout({
      time.tmp <- peakRAM({
        colData(SD)[,1] <- factor(colData(SD)[,1])
        DE.counts <- DESeqDataSetFromMatrix(counts+1,colData(SD),~condition)
        DEseq.results <- DESeq(DE.counts,parallel = T,quiet = T)
      })
      time.mat[i,c("DESeq2","DESeq2_full")] <- time.tmp$Elapsed_Time_sec*12
      total.memory.usage[i,c("DESeq2","DESeq2_full")] <- time.tmp$Total_RAM_Used_MiB
      peak.memory.usage[i,c("DESeq2","DESeq2_full")] <- time.tmp$Peak_RAM_Used_MiB
      
      DEseq.results <- DESeq2::results(DEseq.results,alpha = 0.05)
      adjpvalue.list[["DESeq2"]][i,] <- DEseq.results$padj
      adjpvalue.list[["DESeq2_full"]][i,] <- p.adjust(DEseq.results$pvalue,"BH")
      print("DESeq2 Finish!")
    },timeout = 15000,onTimeout = "error")
  }, error = function(e){
    cat("DESeq2 Error! Simulation ",i,".\n")
    cat("ERROR :",conditionMessage(e),"\n")
  })
  rm(DE.counts,DEseq.results)
  save(adjpvalue.list,time.mat,cluster.acc,cluster.labels.mat,total.memory.usage,peak.memory.usage,file = "./Simulation_V3/NB_200DE_5type_DEG.RData")
  
  ################## EdgeR (filtered) ########################
  library(edgeR)
  tryCatch({
    withTimeout({
      time.tmp <- peakRAM({
        exprSet <- DGEList(counts = as.matrix(counts), group = cluster.labels)
        design <- model.matrix(~cluster.labels)
        colnames(design) <- levels(cluster.labels)
        keep <- filterByExpr(exprSet, design)
        exprSet <- exprSet[keep, , keep.lib.sizes=FALSE]
        exprSet <- calcNormFactors(exprSet)
        exprSet <- estimateDisp(exprSet,design)
        fit <- glmQLFit(exprSet, design)
        qlf <- glmQLFTest(fit, coef=2:(nlevels(cluster.labels)))
      })
      time.mat[i,"EdgeR"] <- time.tmp$Elapsed_Time_sec
      total.memory.usage[i,"EdgeR"] <- time.tmp$Total_RAM_Used_MiB
      peak.memory.usage[i,"EdgeR"] <- time.tmp$Peak_RAM_Used_MiB
      
      colnames(adjpvalue.list[["EdgeR"]]) <- rownames(counts)
      adjpvalue.list[["EdgeR"]][i,rownames(qlf@.Data[[2]])] <- p.adjust(qlf$table$PValue,"BH")
      print("EdgeR Finish!")
    },timeout = 15000,onTimeout = "error")
  }, error = function(e){
    cat("EdgeR filtered Error! Simulation ",i,".\n")
    cat("ERROR :",conditionMessage(e),"\n")
  })
  save(adjpvalue.list,time.mat,cluster.acc,cluster.labels.mat,total.memory.usage,peak.memory.usage,file = "./Simulation_V3/NB_200DE_5type_DEG.RData")
  
  ################## EdgeR (full) ########################
  tryCatch({
    withTimeout({
      time.tmp <- peakRAM({
        exprSet <- DGEList(counts = as.matrix(counts), group = cluster.labels)
        design <- model.matrix(~cluster.labels)
        colnames(design) <- levels(cluster.labels)
        exprSet <- calcNormFactors(exprSet)
        exprSet <- estimateDisp(exprSet,design)
        fit <- glmQLFit(exprSet, design)
        qlf <- glmQLFTest(fit, coef=2:(nlevels(cluster.labels)))
      })
      time.mat[i,"EdgeR_full"] <- time.tmp$Elapsed_Time_sec
      total.memory.usage[i,"EdgeR_full"] <- time.tmp$Total_RAM_Used_MiB
      peak.memory.usage[i,"EdgeR_full"] <- time.tmp$Peak_RAM_Used_MiB
      
      colnames(adjpvalue.list[["EdgeR_full"]]) <- rownames(counts)
      adjpvalue.list[["EdgeR_full"]][i,rownames(qlf@.Data[[2]])] <- p.adjust(qlf$table$PValue,"BH")
      print("EdgeR full Finish!")
    },timeout = 15000,onTimeout = "error")
  }, error = function(e){
    cat("EdgeR full Error! Simulation ",i,".\n")
    cat("ERROR :",conditionMessage(e),"\n")
  })
  save(adjpvalue.list,time.mat,cluster.acc,cluster.labels.mat,total.memory.usage,peak.memory.usage,file = "./Simulation_V3/NB_200DE_5type_DEG.RData")
  rm(exprSet,design,fit,qlf,keep)
  
  ################## MAST ########################
  library(MAST)
  tryCatch({
    withTimeout({
      time.tmp <- peakRAM({
        mast.SD <- FromMatrix(as.matrix(normcounts(SD)),data.frame(condition=colData(SD)[,1],wellKey = colnames(counts)),check_sanity = F)
        zlm.output <- zlm(~ condition,mast.SD)
        zlm.lr <- lrTest(zlm.output, 'condition')
      })
      
      time.mat[i,"MAST"] <- time.tmp$Elapsed_Time_sec
      total.memory.usage[i,"MAST"] <- time.tmp$Total_RAM_Used_MiB
      peak.memory.usage[i,"MAST"] <- time.tmp$Peak_RAM_Used_MiB
      
      MAST.tmp <- zlm.lr[,,'Pr(>Chisq)']
      adjpvalue.list[["MAST"]][i,] <- p.adjust(MAST.tmp[,3],"BH")
      print("MAST Finish!")
    },timeout = 15000,onTimeout = "error")
  }, error = function(e){
    cat("MAST Error! Simulation ",i,".\n")
    cat("ERROR :",conditionMessage(e),"\n")
  })
  save(adjpvalue.list,time.mat,cluster.acc,cluster.labels.mat,total.memory.usage,peak.memory.usage,file = "./Simulation_V3/NB_200DE_5type_DEG.RData")
  rm(mast.SD,zlm.output,zlm.lr,MAST.tmp)
  gc(verbose = F)
  # detach("package:MAST",unload = TRUE)
  
  ################## ZIAQ ########################
  # library(ZIAQ)
  # tryCatch({
  #   withTimeout({
  #     time.tmp <- Sys.time()
  #     ZIAQ.result <- ziaq(normcounts(SD),colDat = colData(SD),log_i = F,parallel = T,no.core = 2)
  #     time.mat[i,"ZIAQ"] <- difftime(Sys.time(),time.tmp,units = "secs")*2
  #     adjpvalue.list[["ZIAQ"]][i,] <- p.adjust(ZIAQ.result[["pvalue"]],"BH")
  #     print("ZIAQ Finish!")
  #   },timeout = 15000,onTimeout = "error")
  # }, error = function(e){
  #   cat("ZIAQ Error! Simulation ",i,".\n")
  #   cat("ERROR :",conditionMessage(e),"\n")
  # })
  # save(adjpvalue.list,time.mat,cluster.acc,cluster.labels.mat,file = "./Simulation_V3/NB_200DE_5type_DEG.RData")
  # # detach("package:ZIAQ",unload = TRUE)
  # rm(ZIAQ.result)
  
  # ################## scDD ########################
  # library(scDD)
  # normcounts(SD) <- as.matrix(normcounts(SD))
  # tryCatch({
  #   withTimeout({
  #     time.tmp <- Sys.time()
  #     SD.tmp <- preprocess(SD, zero.thresh=0.9,scran_norm = T)
  #     prior_param <- list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
  #     scDD.result <- scDD(SD.tmp,prior_param = prior_param,categorize = F,param = snowparam)
  #     scDD.list <- results(scDD.result)
  #     time.mat[i,"scDD"] <- difftime(Sys.time(),time.tmp,units = "secs")*12
  #     adjpvalue.list[["scDD"]][i,] <- scDD.list$nonzero.pvalue.adj
  #     print("scDD Finish!")
  #   },timeout = 15000,onTimeout = "error")
  # }, error = function(e){
  #   cat("scDD Error! Simulation ",i,".\n")
  #   cat("ERROR :",conditionMessage(e),"\n")
  # })
  # save(adjpvalue.list,time.mat,cluster.acc,cluster.labels.mat,file = "./Simulation_V3/NB_200DE_5type_DEG.RData")
  # rm(prior_param,scDD.result,scDD.list,SD.tmp)
  # gc(verbose = F)
  # # detach("package:scDD",unload = TRUE)
  # 
  # ################## sigEMD ########################
  # # SigEMD is based on R package "aod","arm","fdrtool","lars"
  # library(aod)
  # library(arm)
  # library(fdrtool)
  # library(lars)
  # source("./SigEMD-0.21.1/FunImpute.R")
  # source("./SigEMD-0.21.1/SigEMDHur.R")
  # source("./SigEMD-0.21.1/SigEMDnonHur.R")
  # source("./SigEMD-0.21.1/plot_sig.R")
  # tryCatch({
  #   withTimeout({
  #     time.tmp <- Sys.time()
  #     condition.SD <- colData(SD)[,1]
  #     names(condition.SD) <- colnames(normcounts(SD))
  #     results<- calculate_single(data = normcounts(SD),condition =  condition.SD,Hur_gene = NULL, binSize=0.2,nperm=100)
  #     time.mat[i,"sigEMD"] <- difftime(Sys.time(),time.tmp,units = "secs")
  #     emd <- results$emdall
  #     adjpvalue.list[["sigEMD"]][i,] <- emd[,"padjust"]
  #     print("sigEMD Finish!")
  #   },timeout = 15000,onTimeout = "error")
  # }, error = function(e){
  #   cat("sigEMD Error! Simulation ",i,".\n")
  #   cat("ERROR :",conditionMessage(e),"\n")
  # })
  # save(adjpvalue.list,time.mat,cluster.acc,cluster.labels.mat,file = "./Simulation_V3/NB_200DE_5type_DEG.RData")
  # rm(condition.SD,results,emd)
  # gc(verbose = F)
  # 
  
  ################## scDE ########################
  # library(scde)
  # tryCatch({
  #   withTimeout({
  #     time.tmp <- Sys.time()
  #     counts2 <- apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x})
  #     o.ifm2 <- scde.error.models(counts = counts2, groups = colData(SD)[,1], n.cores = 12, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 0)
  #     o.prior <- scde.expression.prior(models = o.ifm2, counts = counts2, length.out = 400, show.plot = FALSE)
  #     ediff <- scde.expression.difference(o.ifm2, counts2, o.prior, groups  =  factor(colData(SD)[,1]), n.randomizations  =  100, n.cores  =  12, verbose  =  0)
  #     time.mat[i,"SCDE"] <- difftime(Sys.time(),time.tmp,units = "secs")*12
  # 
  #     adjpvalue.list[["SCDE"]][i,] <- pnorm(abs(ediff$cZ),lower.tail = F)
  #     print("scDE Finish!")
  #   },timeout = 15000,onTimeout = "error")
  # }, error = function(e){
  #   cat("scDE Error! Simulation ",i,".\n")
  #   cat("ERROR :",conditionMessage(e),"\n")
  # })
  # save(adjpvalue.list,time.mat,cluster.acc,cluster.labels.mat,file = "./Simulation_V3/NB_200DE_5type_DEG.RData")
  # rm(o.ifm2,o.prior,ediff,counts2)
  # detach("package:scde",unload = TRUE)
  # gc(verbose = F)
  
  bpstop(snowparam)
  
  ################## Festem 3groups prior #####################
  for (k in 1:100){
    seurat.SD <- FindClusters(object = seurat.SD, resolution = 0.01*k, verbose = FALSE)
    if (nlevels(seurat.SD@active.ident)>=3){
      print(paste0("Resolution: ", 0.01*k,"; Group number: ",nlevels(seurat.SD@active.ident)))
      break
    }
  }
  cluster.labels <- factor(as.numeric(seurat.SD@active.ident))
  alpha.label <- numeric(nlevels(cluster.labels)-1)
  for (g in 1:length(alpha.label)) {
    alpha.label[g] <- sum(cluster.labels==g)/ncol(seurat.SD)
  }
  cl <- makeCluster(getOption("cl.cores", 12))
  clusterSetRNGStream(cl, iseed = seedlist[i])
  ################## Festem_3group_0.05 ########################
  time.tmp <- peakRAM(
    em.result.3g <- parApply(cl,counts,1,em.stat,
                             alpha.ini=rbind(alpha.label,c(rep(1/nlevels(cluster.labels),length(alpha.label)))),
                             k0=150,C=1e-3,labels = cluster.labels,group.num = nlevels(cluster.labels),
                             prior.weight = 0.05,earlystop = 1e-5)
  )
  time.mat[i,"Festem_3group_0.05"] <- time.tmp$Elapsed_Time_sec*12
  total.memory.usage[i,"Festem_3group_0.05"] <- time.tmp$Total_RAM_Used_MiB
  peak.memory.usage[i,"Festem_3group_0.05"] <- time.tmp$Peak_RAM_Used_MiB
  adjpvalue.list[["Festem_3group_0.05"]][i,] <- p.adjust(em.result.3g[1,],"BH")
  
  ################## Festem_3group_0.05f ########################
  time.tmp <- peakRAM(
    em.result.9.3g <- parApply(cl,counts,1,em.stat,
                               alpha.ini=rbind(alpha.label,c(rep(1/nlevels(cluster.labels),length(alpha.label)))),
                               k0=150,C=1e-3,labels = cluster.labels,group.num = nlevels(cluster.labels),
                               prior.weight = 0.9,earlystop = 1e-5)
  )
  time.mat[i,"Festem_3group_0.05f"] <- time.tmp$Elapsed_Time_sec*12+time.mat[i,"Festem_3group_0.05"]
  total.memory.usage[i,"Festem_3group_0.05f"] <- max(time.tmp$Total_RAM_Used_MiB,total.memory.usage[i,"Festem_3group_0.05"])
  peak.memory.usage[i,"Festem_3group_0.05f"] <- max(time.tmp$Peak_RAM_Used_MiB,peak.memory.usage[i,"Festem_3group_0.05"])
  
  adjpvalue.list[["Festem_3group_0.05f"]][i,] <- adjpvalue.list[["Festem_3group_0.05"]][i,]
  adjpvalue.list[["Festem_3group_0.05f"]][i,em.result.9.3g[2,]<0 & (!is.na(em.result.9.3g[2,])) ] <- NA
  print("Festem Test 3groups Finish!")
  save(adjpvalue.list,time.mat,cluster.acc,cluster.labels.mat,total.memory.usage,peak.memory.usage,file = "./Simulation_V3/NB_200DE_5type_DEG.RData")
  stopCluster(cl)
  gc(verbose = F)
  
  ################## Festem 7groups prior #####################
  for (k in 1:100){
    seurat.SD <- FindClusters(object = seurat.SD, resolution = 0.01*k, verbose = FALSE)
    if (nlevels(seurat.SD@active.ident)>=7){
      print(paste0("Resolution: ", 0.01*k,"; Group number: ",nlevels(seurat.SD@active.ident)))
      break
    }
  }
  cluster.labels <- factor(as.numeric(seurat.SD@active.ident))
  alpha.label <- numeric(nlevels(cluster.labels)-1)
  for (g in 1:length(alpha.label)) {
    alpha.label[g] <- sum(cluster.labels==g)/ncol(seurat.SD)
  }
  cl <- makeCluster(getOption("cl.cores", 12))
  clusterSetRNGStream(cl, iseed = seedlist[i])
  ################## Festem_7group_0.05 ########################
  time.tmp <- peakRAM(
    em.result.4g <- parApply(cl,counts,1,em.stat,
                             alpha.ini=rbind(alpha.label,c(rep(1/nlevels(cluster.labels),length(alpha.label)))),
                             k0=150,C=1e-3,labels = cluster.labels,group.num = nlevels(cluster.labels),
                             prior.weight = 0.05,earlystop = 1e-5)
  )
  time.mat[i,"Festem_7group_0.05"] <- time.tmp$Elapsed_Time_sec*12
  total.memory.usage[i,"Festem_7group_0.05"] <- time.tmp$Total_RAM_Used_MiB
  peak.memory.usage[i,"Festem_7group_0.05"] <- time.tmp$Peak_RAM_Used_MiB
  
  adjpvalue.list[["Festem_7group_0.05"]][i,] <- p.adjust(em.result.4g[1,],"BH")
  
  ################## Festem_7group_0.05f ########################
  time.tmp <- peakRAM(
    em.result.9.4g <- parApply(cl,counts,1,em.stat,
                               alpha.ini=rbind(alpha.label,c(rep(1/nlevels(cluster.labels),length(alpha.label)))),
                               k0=150,C=1e-3,labels = cluster.labels,group.num = nlevels(cluster.labels),
                               prior.weight = 0.9,earlystop = 1e-5)
  )
  time.mat[i,"Festem_7group_0.05f"] <- time.tmp$Elapsed_Time_sec*12+time.mat[i,"Festem_7group_0.05"]
  total.memory.usage[i,"Festem_7group_0.05f"] <- max(time.tmp$Total_RAM_Used_MiB,total.memory.usage[i,"Festem_7group_0.05"])
  peak.memory.usage[i,"Festem_7group_0.05f"] <- max(time.tmp$Peak_RAM_Used_MiB,peak.memory.usage[i,"Festem_7group_0.05"])
  
  adjpvalue.list[["Festem_7group_0.05f"]][i,] <- adjpvalue.list[["Festem_7group_0.05"]][i,]
  adjpvalue.list[["Festem_7group_0.05f"]][i,em.result.9.4g[2,]<0 & (!is.na(em.result.9.4g[2,])) ] <- NA
  print("Festem Test 7groups Finish!")
  save(adjpvalue.list,time.mat,cluster.acc,cluster.labels.mat,total.memory.usage,peak.memory.usage,file = "./Simulation_V3/NB_200DE_5type_DEG.RData")
  stopCluster(cl)
  gc(verbose = F)
  
  ################## Festem 9groups prior #####################
  for (k in 1:100){
    seurat.SD <- FindClusters(object = seurat.SD, resolution = 0.01*k, verbose = FALSE)
    if (nlevels(seurat.SD@active.ident)>=9){
      print(paste0("Resolution: ", 0.01*k,"; Group number: ",nlevels(seurat.SD@active.ident)))
      break
    }
  }
  cluster.labels <- factor(as.numeric(seurat.SD@active.ident))
  alpha.label <- numeric(nlevels(cluster.labels)-1)
  for (g in 1:length(alpha.label)) {
    alpha.label[g] <- sum(cluster.labels==g)/ncol(seurat.SD)
  }
  cl <- makeCluster(getOption("cl.cores", 12))
  clusterSetRNGStream(cl, iseed = seedlist[i])
  ################## Festem_9group_0.05 ########################
  time.tmp <- peakRAM(
    em.result.5g <- parApply(cl,counts,1,em.stat,
                             alpha.ini=rbind(alpha.label,c(rep(1/nlevels(cluster.labels),length(alpha.label)))),
                             k0=150,C=1e-3,labels = cluster.labels,group.num = nlevels(cluster.labels),
                             prior.weight = 0.05,earlystop = 1e-5)
  )
  time.mat[i,"Festem_9group_0.05"] <- time.tmp$Elapsed_Time_sec*12
  total.memory.usage[i,"Festem_9group_0.05"] <- time.tmp$Total_RAM_Used_MiB
  peak.memory.usage[i,"Festem_9group_0.05"] <- time.tmp$Peak_RAM_Used_MiB
  
  adjpvalue.list[["Festem_9group_0.05"]][i,] <- p.adjust(em.result.5g[1,],"BH")
  
  ################## Festem_9group_0.05f ########################
  time.tmp <- peakRAM(
    em.result.9.5g <- parApply(cl,counts,1,em.stat,
                               alpha.ini=rbind(alpha.label,c(rep(1/nlevels(cluster.labels),length(alpha.label)))),
                               k0=150,C=1e-3,labels = cluster.labels,group.num = nlevels(cluster.labels),
                               prior.weight = 0.9,earlystop = 1e-5)
  )
  time.mat[i,"Festem_9group_0.05f"] <- time.tmp$Elapsed_Time_sec*12+time.mat[i,"Festem_9group_0.05"]
  total.memory.usage[i,"Festem_9group_0.05f"] <- max(time.tmp$Total_RAM_Used_MiB,total.memory.usage[i,"Festem_9group_0.05"])
  peak.memory.usage[i,"Festem_9group_0.05f"] <- max(time.tmp$Peak_RAM_Used_MiB,peak.memory.usage[i,"Festem_9group_0.05"])
  
  adjpvalue.list[["Festem_9group_0.05f"]][i,] <- adjpvalue.list[["Festem_9group_0.05"]][i,]
  adjpvalue.list[["Festem_9group_0.05f"]][i,em.result.9.5g[2,]<0 & (!is.na(em.result.9.5g[2,])) ] <- NA
  print("Festem Test 9groups Finish!")
  save(adjpvalue.list,time.mat,cluster.acc,cluster.labels.mat,total.memory.usage,peak.memory.usage,file = "./Simulation_V3/NB_200DE_5type_DEG.RData")
  stopCluster(cl)
  gc(verbose = F)
  
  ################## Feature Selection ###########################
  source("simulation_FS.R")
  save(FS.gene.list,FS.time.mat,file = "./Simulation_V3/NB_200DE_5type_FS.RData")
  gene.list <- list(FS.gene.list[[1]][[i]][1:500],
                    FS.gene.list[[2]][[i]][1:500],
                    FS.gene.list[[3]][[i]][1:500],
                    head(FS.gene.list[[4]][[i]],500),
                    FS.gene.list[[5]][[i]][1:500],
                    FS.gene.list[[6]][[i]][1:500],
                    FS.gene.list[[7]][[i]][1:500],
                    head(FS.gene.list[[8]][[i]],500),
                    head(FS.gene.list[[9]][[i]],500),
                    rownames(seurat.SD),
                    FS.gene.list[[10]][[i]][1:500],
                    FS.gene.list[[11]][[i]][1:500],
                    FS.gene.list[[12]][[i]][1:500])
  umap.list <- vector("list",length(gene.list))
  names(umap.list) <- names(FS.label.list)
  plots.list <- vector("list",length(gene.list))
  truth.tmp <- factor(seurat.SD@meta.data$truth)
  truth.tmp <- as.numeric(truth.tmp)
  
  for (j in 1:length(gene.list)){
    if (length(gene.list[[j]])==0){ next;}
    print(paste0("Method ",j))
    ## kmeans
    clustering <- kmeans(t(seurat.SD@assays$RNA@data[gene.list[[j]],]),2, iter.max = 30, nstart = 5)
    FS.ARI.kmeans.frame[i,j] <- mclust::adjustedRandIndex(clustering$cluster,seurat.SD@meta.data$truth)
    FS.label.kmeans.list[[j]][i,] <- clustering$cluster
    
    ## Louvain
    label.list.tmp <- matrix(nrow = 10, ncol = ncol(seurat.SD))
    ARI.tmp <- numeric(10)
    seurat.SD <- ScaleData(seurat.SD,features = gene.list[[j]])
    seurat.SD <- RunPCA(seurat.SD, verbose = FALSE,features = gene.list[[j]])
    seurat.SD <- FindNeighbors(object = seurat.SD, dims = 1:min(10,length(gene.list[[j]])-1),reduction = "pca")
    
    for (k in 1:10){
      seurat.SD <- FindClusters(object = seurat.SD, resolution = 0.1*k, verbose = FALSE)
      label.list.tmp[k,] <- seurat.SD@active.ident
      num.tmp <- nlevels(seurat.SD@active.ident)
      ARI.tmp[k] <- mclust::adjustedRandIndex(seurat.SD@active.ident,seurat.SD@meta.data$truth)
    }
    
    num.tmp <- abs(num.tmp-5)
    ARI.tmp[num.tmp!=min(num.tmp)] <- 0
    FS.ARI.frame[i,j] <- max(ARI.tmp)
    FS.label.list[[j]][i,] <- label.list.tmp[which.max(ARI.tmp),]
    
    umap.list[[j]] <- RunUMAP(seurat.SD,reduction = "pca", dims = 1:min(10,length(gene.list[[j]]))-1)@reductions[["umap"]]
    dis <- as.matrix(dist(umap.list[[j]]@cell.embeddings, method = "euclidean"))
    FS.SI.frame[i,j] <- summary(cluster::silhouette(x=truth.tmp,dist = dis))[["avg.width"]]
    
    FS.DB.frame[i,j] <- clusterSim::index.DB(umap.list[[j]]@cell.embeddings, 
                                             truth.tmp)$DB
    
    
    umap.tmp <- umap.for.plot(umap.list[[j]],label.list.tmp[which.max(ARI.tmp),])
    
    class_avg <- umap.tmp %>%
      group_by(cluster) %>%
      summarise(
        UMAP_1 = median(UMAP_1),
        UMAP_2 = median(UMAP_2)
      )
    plots.list[[j]] <- ggplot(umap.tmp, aes(x=UMAP_1, y=UMAP_2, color=cluster)) +
      geom_point(cex=0.5) + theme_pubr()+theme(legend.position="none") +
      geom_text(aes(x=UMAP_1,y = UMAP_2,label = cluster), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 4)+
      labs(title = names(umap.list)[j])
    gc(verbose = F)
  }
  save(FS.gene.list,FS.time.mat,FS.ARI.frame,FS.DB.frame,FS.label.list,FS.SI.frame,
       umap.list,plots.list,FS.ARI.kmeans.frame,FS.label.kmeans.list,
       total.memory.usage,peak.memory.usage,file = "./Simulation_V3/NB_200DE_5type_FS.RData")
  gc(verbose = F)
  
  ## True informative gene number
  gene.list <- list(FS.gene.list[[1]][[i]][1:200],
                    FS.gene.list[[2]][[i]][1:200],
                    FS.gene.list[[3]][[i]][1:200],
                    head(FS.gene.list[[4]][[i]],200),
                    FS.gene.list[[5]][[i]][1:200],
                    FS.gene.list[[6]][[i]][1:200],
                    FS.gene.list[[7]][[i]][1:200],
                    head(FS.gene.list[[8]][[i]],200),
                    head(FS.gene.list[[9]][[i]],200),
                    rownames(seurat.SD),
                    FS.gene.list[[10]][[i]][1:200],
                    FS.gene.list[[11]][[i]][1:200],
                    FS.gene.list[[12]][[i]][1:200])
  umap.true.list <- vector("list",length(gene.list))
  names(umap.true.list) <- names(FS.label.true.list)
  plots.true.list <- vector("list",length(gene.list))
  
  for (j in 1:length(gene.list)){
    if (length(gene.list[[j]])==0){ next;}
    
    ## kmeans
    clustering <- kmeans(t(seurat.SD@assays$RNA@data[gene.list[[j]],]),2, iter.max = 30, nstart = 5)
    FS.ARI.kmeans.true.frame[i,j] <- mclust::adjustedRandIndex(clustering$cluster,seurat.SD@meta.data$truth)
    FS.label.kmeans.true.list[[j]][i,] <- clustering$cluster
    
    ## Louvain
    label.list.tmp <- matrix(nrow = 10, ncol = ncol(seurat.SD))
    ARI.tmp <- numeric(10)
    num.tmp <- numeric(10)
    seurat.SD <- ScaleData(seurat.SD,features = gene.list[[j]])
    seurat.SD <- RunPCA(seurat.SD, verbose = FALSE,features = gene.list[[j]])
    seurat.SD <- FindNeighbors(object = seurat.SD, dims = 1:min(10,length(gene.list[[j]])-1),reduction = "pca")
    
    for (k in 1:10){
      seurat.SD <- FindClusters(object = seurat.SD, resolution = 0.1*k, verbose = FALSE)
      label.list.tmp[k,] <- seurat.SD@active.ident
      num.tmp <- nlevels(seurat.SD@active.ident)
      ARI.tmp[k] <- mclust::adjustedRandIndex(seurat.SD@active.ident,seurat.SD@meta.data$truth)
    }
    
    num.tmp <- abs(num.tmp-2)
    ARI.tmp[num.tmp!=min(num.tmp)] <- 0
    FS.ARI.true.frame[i,j] <- max(ARI.tmp)
    FS.label.true.list[[j]][i,] <- label.list.tmp[which.max(ARI.tmp),]
    
    
    umap.true.list[[j]] <- RunUMAP(seurat.SD,reduction = "pca", dims = 1:min(10,length(gene.list[[j]]))-1)@reductions[["umap"]]
    dis <- as.matrix(dist(umap.true.list[[j]]@cell.embeddings, method = "euclidean"))
    FS.SI.true.frame[i,j] <- summary(cluster::silhouette(x=truth.tmp,dist = dis))[["avg.width"]]
    
    FS.DB.true.frame[i,j] <- clusterSim::index.DB(umap.true.list[[j]]@cell.embeddings, 
                                                  truth.tmp)$DB
    
    
    umap.tmp <- umap.for.plot(umap.true.list[[j]],label.list.tmp[which.max(ARI.tmp),])
    
    class_avg <- umap.tmp %>%
      group_by(cluster) %>%
      summarise(
        UMAP_1 = median(UMAP_1),
        UMAP_2 = median(UMAP_2)
      )
    plots.true.list[[j]] <- ggplot(umap.tmp, aes(x=UMAP_1, y=UMAP_2, color=cluster)) +
      geom_point(cex=0.5) + theme_pubr()+theme(legend.position="none") +
      geom_text(aes(x=UMAP_1,y = UMAP_2,label = cluster), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 4)+
      labs(title = names(umap.list)[j])
    gc(verbose = F)
  }
  save(FS.gene.list,FS.time.mat,FS.ARI.frame,FS.DB.frame,FS.label.list,FS.SI.frame,
       umap.list,plots.list,FS.ARI.kmeans.frame,FS.label.kmeans.list,
       FS.ARI.true.frame,FS.DB.true.frame,FS.label.true.list,FS.SI.true.frame,
       umap.true.list,plots.true.list,FS.ARI.kmeans.true.frame,FS.label.kmeans.true.list,
       total.memory.usage,peak.memory.usage,
       file = "./Simulation_V3/NB_200DE_5type_FS.RData")
  gc(verbose = F)
  
  rm(list = ls())
  load("./Simulation_V3/NB_200DE_5type_workspace.RData")
  load("./Simulation_V3/NB_200DE_5type_FS.RData")
  load("./Simulation_V3/NB_200DE_5type_DEG.RData")
}