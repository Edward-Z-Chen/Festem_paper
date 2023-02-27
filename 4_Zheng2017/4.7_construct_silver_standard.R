library(Seurat)
myRunMoransI <- function(data, pos,num.core = 24){
  # Parallel implementation of Seurat's function, RunMoransI.
  # Relying on "Rfast2" package.
  require(parallel)
  cl <- makeCluster(getOption("cl.cores", num.core))
  data <- as.matrix(data)
  pos.dist <- dist(x = pos)
  pos.dist.mat <- as.matrix(x = pos.dist)
  weights <- 1/pos.dist.mat^2
  diag(x = weights) <- 0
  results <- parSapply(cl,X = 1:nrow(x = data), FUN = function(x) {
    Rfast2::moranI(data[x, ], weights)
  })
  pcol <- 2
  results <- data.frame(observed = unlist(x = results[1, ]), 
                        p.value = unlist(x = results[pcol, ]))
  rownames(x = results) <- rownames(x = data)
  stopCluster(cl)
  return(results)
}

# Batch 1 -----------------------------------------------------------------

pbmc <- readRDS("./results/b1.rds")
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc,nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc,features = VariableFeatures(pbmc),verbose = F)
ElbowPlot(pbmc)
pbmc <- RunUMAP(pbmc,dims = 1:30)
umap_hvg <- pbmc@reductions[["umap"]]@cell.embeddings

# Housekeeping genes ------------------------------------------------------
load("Housekeeping_GenesHuman.RData")
Housekeeping_Genes <- dplyr::filter(Housekeeping_Genes,Gene.name%in%rownames(pbmc))
genes <- Housekeeping_Genes$Gene.name
genes <- unique(genes)
moran_h <- myRunMoransI(pbmc@assays$RNA@data[genes,],umap_hvg)

# Non-housekeeping genes --------------------------------------------------
genes <- setdiff(rownames(pbmc),Housekeeping_Genes$Gene.name)
moran_nonh <- myRunMoransI(pbmc@assays$RNA@data[genes,],umap_hvg)


# Non-zero expression percentage ------------------------------------------
cal_express_percent <- function(x){
  sum(x>0)/length(x)
}
gene_percent <- apply(pbmc@assays$RNA@counts,1,cal_express_percent)

save(moran_h,moran_nonh,gene_percent,file = "./results/b1_silver_standard.RData")

# Batch 2 -----------------------------------------------------------------

pbmc <- readRDS("./results/b2.rds")
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc,nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc,features = VariableFeatures(pbmc),verbose = F)
ElbowPlot(pbmc)
pbmc <- RunUMAP(pbmc,dims = 1:30)
umap_hvg <- pbmc@reductions[["umap"]]@cell.embeddings

# Housekeeping genes ------------------------------------------------------
load("Housekeeping_GenesHuman.RData")
Housekeeping_Genes <- dplyr::filter(Housekeeping_Genes,Gene.name%in%rownames(pbmc))
genes <- Housekeeping_Genes$Gene.name
genes <- unique(genes)
moran_h <- myRunMoransI(pbmc@assays$RNA@data[genes,],umap_hvg)

# Non-housekeeping genes --------------------------------------------------
genes <- setdiff(rownames(pbmc),Housekeeping_Genes$Gene.name)
moran_nonh <- myRunMoransI(pbmc@assays$RNA@data[genes,],umap_hvg)


# Non-zero expression percentage ------------------------------------------
cal_express_percent <- function(x){
  sum(x>0)/length(x)
}
gene_percent <- apply(pbmc@assays$RNA@counts,1,cal_express_percent)

save(moran_h,moran_nonh,gene_percent,file = "./results/b2_silver_standard.RData")
