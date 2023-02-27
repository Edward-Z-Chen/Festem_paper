library(ggplot2)
library(tidyverse)
library(scales)
library(Seurat)
library(harmony)
library(scattermore)
nonepi <- readRDS("NonEpi.rds")
load("./results/iCCA_hvggenes.RData")
nonepi <- NormalizeData(nonepi)
umap.for.plot <- function(umap,cluster){
  umap <- umap@cell.embeddings
  umap <- as.data.frame(umap)
  cbind(umap,cluster = cluster)
}


gene.list <- list(EM[1:7000],
                  hvgvst[1:7000],
                  hvgdisp[1:7000],
                  dub,
                  devianceFS[1:7000],
                  trendvar[1:7000])
umap.list <- vector("list",length(gene.list))
label.list <- vector("list",length(gene.list))
names(label.list) <- c("Festem","HVGvst","HVGdisp","DUBStepR","devianceFS","trendVar")
names(umap.list) <- c("Festem","HVGvst","HVGdisp","DUBStepR","devianceFS","trendVar")
plots.list <- vector("list",length(gene.list))
for (i in 1:length(gene.list)){
  nonepi <- ScaleData(nonepi,features = gene.list[[i]])
  nonepi <- RunPCA(nonepi, verbose = FALSE,features = gene.list[[i]])
  nonepi <- RunHarmony(nonepi,"Patient",plot_convergence = T)
  nonepi <- FindNeighbors(object = nonepi, dims = 1:30,reduction = "harmony")
  nonepi <- FindClusters(object = nonepi, resolution = 0.6)
  label.list[[i]] <- nonepi@active.ident
  umap.list[[i]] <- RunUMAP(nonepi,reduction = "harmony", dims = 1:30)@reductions[["umap"]]
  umap.tmp <- umap.for.plot(umap.list[[i]],label.list[[i]])
  
  class_avg <- umap.tmp %>%
    group_by(cluster) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  plots.list[[i]] <- ggplot(umap.tmp, aes(x=UMAP_1, y=UMAP_2, color=cluster)) + 
    geom_scattermore(pointsize = 1)+
    theme_bw()+theme(legend.position="none") +
    geom_text(aes(x=UMAP_1,y = UMAP_2,label = cluster), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 4)+
    labs(title = names(label.list)[i])
}

save(gene.list,label.list,umap.list,plots.list,file = "./results/iCCA_clustering_UMAP.RData")


# # Different resolution ----------------------------------------------------
# ref_label <- label.list[[1]]
# label.list <- vector("list",length(gene.list))
# names(label.list) <- c("Festem","HVGvst","HVGDisp","DUBStepR","devianceFS","trendVar")
# label_ARI <- matrix(nrow = length(gene.list),ncol = 15)
# total_num <- matrix(nrow = length(gene.list),ncol = 15)
# rownames(label_ARI) <- names(label.list)
# for (i in 1:length(gene.list)){
#   label.list[[i]] <- matrix(nrow = 15, ncol = ncol(nonepi))
#   nonepi <- ScaleData(nonepi,features = gene.list[[i]])
#   nonepi <- RunPCA(nonepi, verbose = FALSE,features = gene.list[[i]])
#   nonepi <- FindNeighbors(object = nonepi, dims = 1:15)
#   for (j in 1:15){
#     nonepi <- FindClusters(object = nonepi, resolution = 0.1*j)
#     label.list[[i]][j,] <- nonepi@active.ident
#     label_ARI[i,j] <- mclust::adjustedRandIndex(nonepi@active.ident,ref_label)
#     total_num[i,j] <- nlevels(nonepi@active.ident)
#   }
# }
# save(label.list,label_ARI,total_num,file = "./results/nonepi3k_resolution_ARI.RData")
