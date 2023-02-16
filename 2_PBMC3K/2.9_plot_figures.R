library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(scales)
library(RColorBrewer)
my.color <- hue_pal()(13)
names(my.color) <- 1:13
umap.for.plot <- function(umap,cluster){
  umap <- umap@cell.embeddings
  umap <- as.data.frame(umap)
  cbind(umap,cluster = cluster)
}

# Figure 3 (A) -- left top ------------------------------------------------
## Summarising results from various DEG detection methods
results <- matrix(nrow = nrow(pbmc),ncol = 13,
                  dimnames = list(rownames(pbmc),
                                  c("Festem","DEseq2", "DEseq2-full",       
                                    "EdgeR", "EdgeR-full","MAST",              
                                    "Wilcoxon","MAST-f",            
                                    "Wilcoxon-f","FC",
                                    "singleCellHaystack-PCA",
                                    "singleCellHaystack-UMAP","TN_test")))
pbmc <- readRDS("./results/pbmc3k.rds")
pbmc <- as.SingleCellExperiment(pbmc)
load("./results/pbmc3k_label.RData")
cluster.label[c("CGGGCATGACCCAA-1","CTTGATTGATCTTC-1")] <- "Platelet"
pbmc <- pbmc[,cluster.label!="Platelet"]
cluster.label <- factor(cluster.label[cluster.label!="Platelet"])
pbmc <- AddMetaData(pbmc,cluster.label,"HVG")
FC.list <- vector("list",8)
for (i in 1:8){
  FC.list[[i]] <- FoldChange(pbmc,ident.1 = levels(cluster.label)[i],
                             group.by = "HVG")[,1]
}
FC.list <- matrix(unlist(FC.list),ncol = 8,byrow = F)
FC <- apply(FC.list,1,function(x){max(abs(x))})
results[,"FC"] <- FC
load("./results/pbmc_Festem.RData")
results[colnames(em.result),"Festem"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result9)[em.result9[2,]<0],"Festem"] <- NA
load("./results/pbmc_DEseq2.RData")
results[DEseq.results@rownames,"DESeq2"] <- DEseq.results@listData$padj
results[DEseq.results@rownames,"DESeq2-full"] <- p.adjust(DEseq.results@listData$pvalue,"BH")
load("./results/pbmc_edgeR.RData")
results[names(EdgeR.result),"EdgeR"] <- p.adjust(EdgeR.result,"BH")
load("./results/pbmc_edgeR_full.RData")
results[names(EdgeR.result),"EdgeR-full"] <- p.adjust(EdgeR.result,"BH")
load("./results/pbmc_mast.RData")
results[rownames(mast.results),"MAST"] <- p.adjust(mast.results[,3],"BH")
load("./results/pbmc_wilcoxon.RData")
wil.result <- matrix(unlist(wil.result),ncol = 8,byrow = F)
wil.result <- apply(wil.result,1,function(x){min(x)*8})
results[,"Wilcoxon"] <- p.adjust(wil.result,"BH")
results[,"Wilcoxon-f"] <- results[,"Wilcoxon"]
results[results[,"FC"]<=0.2,"Wilcoxon-f"] <- NA
results[,"MAST-f"] <- results[,"MAST"]
results[results[,"FC"]<=0.2,"MAST-f"] <- NA
load("./results/pbmc_haystack.RData")
results[,"singleCellHaystack-PCA"] <- exp(haystack_pca$results$log.p.adj)
results[,"singleCellHaystack-UMAP"] <- exp(haystack_umap$results$log.p.adj)

tn_test <- read.csv("./results/pbmc_TN_test.csv",
                    header = F)
tn_test <- as.matrix(tn_test)
tn_test <- apply(tn_test,1,function(x){min(x)*length(x)})
results[,"TN_test"] <- p.adjust(tn_test,"BH")
save(results,file = "./results/pbmc3k_DEG_results.RData")
## Plots