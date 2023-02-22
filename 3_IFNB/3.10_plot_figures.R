library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(scales)
library(RColorBrewer)
my.color <- hue_pal()(17)
names(my.color) <- 1:17
tsne.for.plot <- function(tsne,cluster){
  tsne <- tsne@cell.embeddings
  tsne <- as.data.frame(tsne)
  cbind(tsne,cluster = cluster)
}
if (!file.exists("./figures")){
  dir.create("./figures")
}

# Figure 3 (A) -- top right and Figure S5 second from top ------------------------------------------------
## Summarizing results from various DEG detection methods
ifnb <- readRDS("./results/ifnb_ctrl.rds")
results <- matrix(nrow = nrow(ifnb),ncol = 19,
                  dimnames = list(rownames(ifnb),
                                  c("Festem","Festem_gamma0.01","Festem_13g","Festem_17g",
                                    "DEseq2", "DEseq2-full",       
                                    "EdgeR", "EdgeR-full","MAST",              
                                    "Wilcoxon","MAST-f",            
                                    "Wilcoxon-f","FC",
                                    "singleCellHaystack-PCA",
                                    "singleCellHaystack-TSNE",
                                    "singleCellHaystack-PCA-20pc",
                                    "singleCellHaystack-TSNE-20pc","TN_test","ROSeq")))
FC.list <- vector("list",nlevels(ifnb@meta.data$seurat_annotations))
for (i in 1:length(FC.list)){
  FC.list[[i]] <- FoldChange(ifnb,ident.1 = levels(ifnb@meta.data$seurat_annotations)[i],
                             group.by = "seurat_annotations")[,1]
}
FC.list <- matrix(unlist(FC.list),ncol = nlevels(ifnb@meta.data$seurat_annotations),byrow = F)
FC <- apply(FC.list,1,function(x){max(abs(x))})
results[,"FC"] <- FC
load("./results/ifnb_ctrl_Festem.RData")
results[colnames(em.result),"Festem"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result.9)[em.result.9[2,]<0],"Festem"] <- NA
load("./results/ifnb_ctrl_Festem_gamma0.01.RData")
results[colnames(em.result),"Festem_gamma0.01"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result.9)[em.result.9[2,]<0],"Festem_gamma0.01"] <- NA
load("./results/ifnb_ctrl_Festem_13g.RData")
results[colnames(em.result),"Festem_13g"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result.9)[em.result.9[2,]<0],"Festem_13g"] <- NA
load("./results/ifnb_ctrl_Festem_17g.RData")
results[colnames(em.result),"Festem_17g"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result.9)[em.result.9[2,]<0],"Festem_17g"] <- NA
load("./results/ifnb_DEseq2.RData")
results[DEseq.results@rownames,"DEseq2"] <- DEseq.results@listData$padj
results[DEseq.results@rownames,"DEseq2-full"] <- p.adjust(DEseq.results@listData$pvalue,"BH")
load("./results/ifnb_edgeR.RData")
results[names(EdgeR.result),"EdgeR"] <- p.adjust(EdgeR.result,"BH")
load("./results/ifnb_edgeR_full.RData")
results[names(EdgeR.result),"EdgeR-full"] <- p.adjust(EdgeR.result,"BH")
load("./results/ifnb_mast.RData")
results[rownames(mast.results),"MAST"] <- p.adjust(mast.results[,3],"BH")
load("./results/ifnb_wilcoxon.RData")
wil.result <- matrix(unlist(wil.result),ncol = length(wil.result),byrow = F)
wil.result <- apply(wil.result,1,function(x){min(x)*length(wil.result)})
results[,"Wilcoxon"] <- p.adjust(wil.result,"BH")
results[,"Wilcoxon-f"] <- results[,"Wilcoxon"]
results[results[,"FC"]<=0.2,"Wilcoxon-f"] <- NA
results[,"MAST-f"] <- results[,"MAST"]
results[results[,"FC"]<=0.2,"MAST-f"] <- NA
load("./results/ifnb_haystack.RData")
results[,"singleCellHaystack-PCA"] <- exp(haystack_pca$results$log.p.adj)
results[,"singleCellHaystack-TSNE"] <- exp(haystack_umap$results$log.p.adj)
load("./results/ifnb_haystack_20PC.RData")
results[,"singleCellHaystack-PCA-20pc"] <- exp(haystack_pca$results$log.p.adj)
results[,"singleCellHaystack-TSNE-20pc"] <- exp(haystack_umap$results$log.p.adj)
load("./results/ifnb_ROSeq.RData")
results[,"ROSeq"] <- p.adjust(roseq.tmp,"BH")

tn_test <- read.csv("./results/ifnb_TN_test.csv",
                    header = F)
tn_test <- as.matrix(tn_test)
tn_test <- apply(tn_test,1,function(x){min(x)*length(x)})
results[,"TN_test"] <- p.adjust(tn_test,"BH")
save(results,file = "./results/ifnb_ctrl_DEG_results.RData")
