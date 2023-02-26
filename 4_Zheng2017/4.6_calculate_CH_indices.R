library(Seurat)
library(cluster)
library(fpc)

# Batch 1 -----------------------------------------------------------------
pbmc <- readRDS("./results/b1.rds")
load("./results/b1_hvggenes.RData")
pbmc <- NormalizeData(pbmc)
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]


gene.list <- list(EM[1:3000],
                  hvgvst[1:3000],
                  hvgdisp[1:3000],
                  dub,
                  devianceFS[1:3000],
                  trendvar[1:3000])
CH.index <- rep(NA,length(gene.list))

for (j in 1:length(CH.index)){
  if (length(gene.list[[j]])>=3000){
    gene.tmp <- gene.list[[j]]
    pbmc <- ScaleData(pbmc,features = gene.tmp)
    pbmc <- RunPCA(pbmc, verbose = FALSE,features = gene.tmp)
    pbmc <- FindNeighbors(object = pbmc, dims = 1:20)
    pbmc <- FindClusters(object = pbmc, resolution = 0.6)
    pca <- pbmc@reductions$pca@cell.embeddings[,1:20]
    CH.index[j] <- calinhara(x=pca,clustering = as.numeric(pbmc@active.ident))
  }
}
print(CH.index)
save(CH.index,file = "./results/b1_CH.RData")

# Batch 2 -----------------------------------------------------------------
pbmc <- readRDS("./results/b2.rds")
load("./results/b2_hvggenes.RData")
pbmc <- NormalizeData(pbmc)
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]


gene.list <- list(EM[1:3000],
                  hvgvst[1:3000],
                  hvgdisp[1:3000],
                  dub,
                  devianceFS[1:3000],
                  trendvar[1:3000])
CH.index <- rep(NA,length(gene.list))

for (j in 1:length(CH.index)){
  if (length(gene.list[[j]])>=3000){
    gene.tmp <- gene.list[[j]]
    pbmc <- ScaleData(pbmc,features = gene.tmp)
    pbmc <- RunPCA(pbmc, verbose = FALSE,features = gene.tmp)
    pbmc <- FindNeighbors(object = pbmc, dims = 1:20)
    pbmc <- FindClusters(object = pbmc, resolution = 0.6)
    pca <- pbmc@reductions$pca@cell.embeddings[,1:20]
    CH.index[j] <- calinhara(x=pca,clustering = as.numeric(pbmc@active.ident))
  }
}
print(CH.index)
save(CH.index,file = "./results/b2_CH.RData")