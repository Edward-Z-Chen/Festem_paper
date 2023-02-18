# Download preprocessed data -----------------------------------------------------------
library(Seurat)
library(SeuratData)
data("ifnb")
ifnb <- subset(ifnb,subset = orig.ident=="IMMUNE_CTRL")
saveRDS(ifnb,file = "./results/ifnb_ctrl.rds")


# Preclustering -----------------------------------------------------------
ifnb <- readRDS("./results/ifnb_ctrl.rds")
ifnb <- ScaleData(ifnb, features = rownames(ifnb))
ifnb <- RunPCA(ifnb, features = rownames(ifnb),verbose = F)
ifnb <- FindNeighbors(ifnb, dims = 1:20)
ifnb <- FindClusters(ifnb, resolution = 0.9)
ifnb <- RunTSNE(ifnb,dims = 1:20)
DimPlot(ifnb, reduction = "tsne",label = T)

cluster.labels <- ifnb@meta.data$seurat_clusters
levels(cluster.labels) <- 1:nlevels(cluster.labels)
save(cluster.labels,file = "./results/ifnb_ctrl_preclustering.RData")

ifnb <- FindClusters(ifnb, resolution = 0.7)
cluster.labels <- ifnb@meta.data$seurat_clusters
levels(cluster.labels) <- 1:nlevels(cluster.labels)
save(cluster.labels,file = "./results/ifnb_ctrl_preclustering_13g.RData")

ifnb <- FindClusters(ifnb, resolution = 1.2)
cluster.labels <- ifnb@meta.data$seurat_clusters
levels(cluster.labels) <- 1:nlevels(cluster.labels)
save(cluster.labels,file = "./results/ifnb_ctrl_preclustering_17g.RData")