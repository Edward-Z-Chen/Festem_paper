library(Seurat)

# Batch 1 -----------------------------------------------------------------

pbmc <- readRDS("./results/b1.rds")
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc,nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc,features = VariableFeatures(pbmc),verbose = F)
ElbowPlot(pbmc)
pbmc <- RunUMAP(pbmc,dims = 1:20)
umap_hvg <- pbmc@reductions[["umap"]]@cell.embeddings

# Housekeeping genes ------------------------------------------------------
load("Housekeeping_GenesHuman.RData")
Housekeeping_Genes <- dplyr::filter(Housekeeping_Genes,Gene.name%in%rownames(pbmc))
genes <- Housekeeping_Genes$Gene.name
genes <- unique(genes)
moran_h <- RunMoransI(pbmc@assays$RNA@data[genes,],umap_hvg)

# Non-housekeeping genes --------------------------------------------------
genes <- intersect(rownames(pbmc),Housekeeping_Genes$Gene.name)
moran_nonh <- RunMoransI(pbmc@assays$RNA@data[genes,],umap_hvg)


# Non-zero expression percentage ------------------------------------------
cal_express_percent <- function(x){
  sum(x>0)/length(x)
}
gene_percent <- apply(pbmc@assays$RNA@counts,1,cal_express_percent)

save(moran_h,moran_nonh,gene_percent,file = "./results/b2_silver_standard.RData")

# Batch 2 -----------------------------------------------------------------

pbmc <- readRDS("./results/b2.rds")
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc,nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc,features = VariableFeatures(pbmc),verbose = F)
ElbowPlot(pbmc)
pbmc <- RunUMAP(pbmc,dims = 1:20)
umap_hvg <- pbmc@reductions[["umap"]]@cell.embeddings

# Housekeeping genes ------------------------------------------------------
load("Housekeeping_GenesHuman.RData")
Housekeeping_Genes <- dplyr::filter(Housekeeping_Genes,Gene.name%in%rownames(pbmc))
genes <- Housekeeping_Genes$Gene.name
genes <- unique(genes)
moran_h <- RunMoransI(pbmc@assays$RNA@data[genes,],umap_hvg)

# Non-housekeeping genes --------------------------------------------------
genes <- intersect(rownames(pbmc),Housekeeping_Genes$Gene.name)
moran_nonh <- RunMoransI(pbmc@assays$RNA@data[genes,],umap_hvg)


# Non-zero expression percentage ------------------------------------------
cal_express_percent <- function(x){
  sum(x>0)/length(x)
}
gene_percent <- apply(pbmc@assays$RNA@counts,1,cal_express_percent)

save(moran_h,moran_nonh,gene_percent,file = "./results/b2_silver_standard.RData")
