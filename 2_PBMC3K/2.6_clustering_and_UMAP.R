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
load("./results/pbmc3k_hvggenes.RData")
pbmc <- readRDS("./results/pbmc3k.rds")
ref <- pbmc@active.ident
ref[c("CGGGCATGACCCAA-1","CTTGATTGATCTTC-1")] <- "Platelet"
pbmc <- pbmc[,ref!="Platelet"]
gene.list <- list(EM[1:1000],
                  hvgvst[1:1000],
                  hvgdisp[1:1000],
                  dub,
                  devianceFS[1:1000],
                  trendvar[1:1000])
umap.list <- vector("list",length(gene.list))
label.list <- vector("list",length(gene.list))
names(label.list) <- c("Festem","HVGvst","HVGdisp","DUBStepR","devianceFS","trendVar")
names(umap.list) <- c("Festem","HVGvst","HVGdisp","DUBStepR","devianceFS","trendVar")
plots.list <- vector("list",length(gene.list))
for (i in 1:length(gene.list)){
  pbmc <- ScaleData(pbmc,features = gene.list[[i]])
  pbmc <- RunPCA(pbmc, verbose = FALSE,features = gene.list[[i]])
  pbmc <- FindNeighbors(object = pbmc, dims = 1:15)
  pbmc <- FindClusters(object = pbmc, resolution = 1)
  label.list[[i]] <- pbmc@active.ident
  umap.list[[i]] <- RunUMAP(pbmc, reduction = "pca", dims = 1:15)@reductions[["umap"]]
  umap.tmp <- umap.for.plot(umap.list[[i]],label.list[[i]])
  
  class_avg <- umap.tmp %>%
    group_by(cluster) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  plots.list[[i]] <- ggplot(umap.tmp, aes(x=UMAP_1, y=UMAP_2, color=cluster)) + 
    geom_point(cex=0.5) + theme_bw()+theme(legend.position="none") +
    geom_text(aes(x=UMAP_1,y = UMAP_2,label = cluster), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 4)+
    labs(title = names(label.list)[i])
}

save(gene.list,label.list,umap.list,plots.list,file = "./results/pbmc3k_clustering_UMAP.RData")


# Different resolution ----------------------------------------------------
ref_label <- label.list[[1]]
label.list <- vector("list",length(gene.list))
names(label.list) <- c("Festem","HVGvst","HVGDisp","DUBStepR","devianceFS","trendVar")
label_ARI <- matrix(nrow = length(gene.list),ncol = 15)
total_num <- matrix(nrow = length(gene.list),ncol = 15)
rownames(label_ARI) <- names(label.list)
for (i in 1:length(gene.list)){
  label.list[[i]] <- matrix(nrow = 15, ncol = ncol(pbmc))
  pbmc <- ScaleData(pbmc,features = gene.list[[i]])
  pbmc <- RunPCA(pbmc, verbose = FALSE,features = gene.list[[i]])
  pbmc <- FindNeighbors(object = pbmc, dims = 1:15)
  for (j in 1:15){
    pbmc <- FindClusters(object = pbmc, resolution = 0.1*j)
    label.list[[i]][j,] <- pbmc@active.ident
    label_ARI[i,j] <- mclust::adjustedRandIndex(pbmc@active.ident,ref_label)
    total_num[i,j] <- nlevels(pbmc@active.ident)
  }
}
save(label.list,label_ARI,total_num,file = "./results/pbmc3k_resolution_ARI.RData")
