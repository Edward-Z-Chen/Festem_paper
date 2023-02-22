# Download preprocessed data -----------------------------------------------------------
library(Seurat)
library(SeuratData)
InstallData("ifnb")
data("ifnb")
ifnb <- subset(ifnb,subset = orig.ident=="IMMUNE_CTRL")
saveRDS(ifnb,file = "./results/ifnb_ctrl.rds")


# Preclustering -----------------------------------------------------------
ifnb <- readRDS("./results/ifnb_ctrl.rds")
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb,nfeatures = 8000)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb,verbose = F)
ifnb <- FindNeighbors(ifnb, dims = 1:20)
ifnb <- FindClusters(ifnb, resolution = 0.7)
ifnb <- RunTSNE(ifnb,dims = 1:20)
DimPlot(ifnb, reduction = "tsne",label = T)

cluster.labels <- ifnb@meta.data$seurat_clusters
levels(cluster.labels) <- 1:nlevels(cluster.labels)
save(cluster.labels,file = "./results/ifnb_ctrl_preclustering.RData")

ifnb <- FindClusters(ifnb, resolution = 0.6)
cluster.labels <- ifnb@meta.data$seurat_clusters
levels(cluster.labels) <- 1:nlevels(cluster.labels)
save(cluster.labels,file = "./results/ifnb_ctrl_preclustering_13g.RData")

ifnb <- FindClusters(ifnb, resolution = 1)
cluster.labels <- ifnb@meta.data$seurat_clusters
levels(cluster.labels) <- 1:nlevels(cluster.labels)
save(cluster.labels,file = "./results/ifnb_ctrl_preclustering_17g.RData")

# Count matrix for TN test (python) ---------------------------------------
counts <- as.matrix(pbmc@assays$RNA@counts)
save(counts,file = "./results/ifnb_counts_forpython.RData")