library(Seurat)
library(cluster)
library(fpc)
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
CH.index <- rep(NA,length(gene.list))

for (j in 1:length(CH.index)){
  if (length(gene.list[[j]])>=2500){
    gene.tmp <- gene.list[[j]]
    ifnb <- ScaleData(ifnb,features = gene.tmp)
    ifnb <- RunPCA(ifnb, verbose = FALSE,features = gene.tmp)
    ifnb <- FindNeighbors(object = ifnb, dims = 1:25)
    ifnb <- FindClusters(object = ifnb, resolution = 1.5)
    pca <- ifnb@reductions$pca@cell.embeddings[,1:25]
    CH.index[j] <- calinhara(x=pca,clustering = as.numeric(ifnb@active.ident))
  }
}
save(CH.index,file = "./results/ifnb_CH.RData")