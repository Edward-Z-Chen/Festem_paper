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
load("./results/ifnb_ctrl_hvggenes.RData")
ifnb <- readRDS("./results/ifnb_ctrl.rds")
ifnb <- NormalizeData(ifnb)
ifnb <- ifnb[,!ifnb@meta.data$seurat_annotations%in%c("Eryth","Mk")]
gene.list <- list(EM[1:2500],
                  hvgvst[1:2500],
                  hvgdisp[1:2500],
                  dub,
                  devianceFS[1:2500],
                  trendvar[1:2500])
tsne.list <- vector("list",length(gene.list))
label.list <- vector("list",length(gene.list))
names(label.list) <- c("Festem","HVGvst","HVGdisp","DUBStepR","devianceFS","trendVar")
names(tsne.list) <- c("Festem","HVGvst","HVGdisp","DUBStepR","devianceFS","trendVar")
plots.list <- vector("list",length(gene.list))
for (i in 1:length(gene.list)){
  ifnb <- ScaleData(ifnb,features = gene.list[[i]])
  ifnb <- RunPCA(ifnb, verbose = FALSE,features = gene.list[[i]])
  ifnb <- FindNeighbors(object = ifnb, dims = 1:min(25,length(gene.list[[i]])-1))
  ifnb <- FindClusters(object = ifnb, resolution = 1.5)
  label.list[[i]] <- ifnb@active.ident
  tsne.list[[i]] <- RunTSNE(ifnb, reduction = "pca", dims = 1:min(25,length(gene.list[[i]])-1))@reductions[["tsne"]]
  tsne.tmp <- tsne.for.plot(tsne.list[[i]],label.list[[i]])
  
  class_avg <- tsne.tmp %>%
    group_by(cluster) %>%
    summarise(
      tSNE_1 = median(tSNE_1),
      tSNE_2 = median(tSNE_2)
    )
  plots.list[[i]] <- ggplot(tsne.tmp, aes(x=tSNE_1, y=tSNE_2, color=cluster)) + 
    geom_point(cex=0.5) + theme_bw()+theme(legend.position="none") +
    geom_text(aes(x=tSNE_1,y = tSNE_2,label = cluster), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 4)+
    labs(title = names(label.list)[i])
}

save(gene.list,label.list,tsne.list,plots.list,file = "./results/ifnb_ctrl_clustering_tSNE.RData")


# Different resolution ----------------------------------------------------
ref_label <- label.list[[1]]
label.list <- vector("list",length(gene.list))
names(label.list) <- c("Festem","HVGvst","HVGDisp","DUBStepR","devianceFS","trendVar")
label_ARI <- matrix(nrow = length(gene.list),ncol = 15)
total_num <- matrix(nrow = length(gene.list),ncol = 15)
rownames(label_ARI) <- names(label.list)
for (i in 1:length(gene.list)){
  label.list[[i]] <- matrix(nrow = 15, ncol = ncol(ifnb))
  ifnb <- ScaleData(ifnb,features = gene.list[[i]])
  ifnb <- RunPCA(ifnb, verbose = FALSE,features = gene.list[[i]])
  ifnb <- FindNeighbors(object = ifnb, dims = 1:min(25,length(gene.list[[i]])-1))
  for (j in 1:15){
    ifnb <- FindClusters(object = ifnb, resolution = 0.1*j+0.5)
    label.list[[i]][j,] <- ifnb@active.ident
    label_ARI[i,j] <- mclust::adjustedRandIndex(ifnb@active.ident,ref_label)
    total_num[i,j] <- nlevels(ifnb@active.ident)
  }
}
save(label.list,label_ARI,total_num,file = "./results/ifnb_ctrl_resolution_ARI.RData")
