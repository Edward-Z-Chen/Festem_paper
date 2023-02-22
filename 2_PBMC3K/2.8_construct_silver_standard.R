library(Seurat)
pbmc <- readRDS("./results/pbmc3k.rds")
load("./results/pbmc3k_clustering_UMAP.RData")
pbmc <- pbmc[,rownames(plots.list[[1]]$data)]
pbmc <- NormalizeData(pbmc)
umap_hvg <- umap.list[[2]]@cell.embeddings
umap_hvg <- umap_hvg[rownames(plots.list[[2]]$data),]

# Housekeeping genes ------------------------------------------------------
load("Housekeeping_GenesHuman.RData")
Housekeeping_Genes <- dplyr::filter(Housekeeping_Genes,Gene.name%in%rownames(pbmc))
genes <- Housekeeping_Genes$Gene.name
moran_h <- RunMoransI(pbmc@assays$RNA@data[genes,],umap_hvg)

# Non-housekeeping genes --------------------------------------------------
genes <- intersect(rownames(pbmc),Housekeeping_Genes$Gene.name)
moran_nonh <- RunMoransI(pbmc@assays$RNA@data[genes,],umap_hvg)


# Non-zero expression percentage ------------------------------------------
cal_express_percent <- function(x){
  sum(x>0)/length(x)
}
gene_percent <- apply(pbmc@assays$RNA@counts,1,cal_express_percent)

save(moran_h,moran_nonh,gene_percent,file = "./results/pbmc3k_silver_standard.RData")
