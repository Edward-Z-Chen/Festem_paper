library(Seurat)
ifnb <- readRDS("./results/ifnb_ctrl.rds")
load("./results/ifnb_ctrl_clustering_tSNE.RData")

ifnb <- ifnb[,rownames(plots.list[[2]]$data)]
ifnb <- NormalizeData(ifnb)
umap_hvg <- umap.list[[2]]@cell.embeddings
umap_hvg <- umap_hvg[rownames(plots.list[[2]]$data),]

# Housekeeping genes ------------------------------------------------------
load("Housekeeping_GenesHuman.RData")
Housekeeping_Genes <- dplyr::filter(Housekeeping_Genes,Gene.name%in%rownames(ifnb))
moran_h <- RunMoransI(ifnb@assays$RNA@data[genes,],umap_hvg)

# Non-housekeeping genes --------------------------------------------------
genes <- intersect(rownames(ifnb),Housekeeping_Genes$Gene.name)
moran_nonh <- RunMoransI(ifnb@assays$RNA@data[genes,],umap_hvg)


# Non-zero expression percentage ------------------------------------------
cal_express_percent <- function(x){
  sum(x>0)/length(x)
}
gene_percent <- apply(ifnb@assays$RNA@counts,1,cal_express_percent)

save(moran_h,moran_nonh,gene_percent,file = "./results/ifnb_ctrl_silver_standard.RData")
