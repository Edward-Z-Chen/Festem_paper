library(Seurat)
library(cluster)
library(fpc)
load("./results/pbmc3k_hvggenes.RData")
pbmc <- readRDS("./results/pbmc3k.rds")
ref <- pbmc@active.ident
ref[c("CGGGCATGACCCAA-1","CTTGATTGATCTTC-1")] <- "Platelet"
pbmc <- pbmc[,ref!="Platelet"]

gene.list <- list(EM[1:1000],hvgvst[1:1000],hvgdisp[1:1000],dub,
                  devianceFS[1:1000],trendvar[1:1000])
CH.index <- rep(NA,length(gene.list))

for (j in 1:length(CH.index)){
  if (length(gene.list[[j]])>=1000){
    gene.tmp <- gene.list[[j]]
    pbmc <- ScaleData(pbmc,features = gene.tmp)
    pbmc <- RunPCA(pbmc, verbose = FALSE,features = gene.tmp)
    pbmc <- FindNeighbors(object = pbmc, dims = 1:15)
    pbmc <- FindClusters(object = pbmc, resolution = 1)
    pca <- pbmc@reductions$pca@cell.embeddings[,1:15]
    CH.index[j] <- calinhara(x=pca,clustering = as.numeric(label.list[[i]]))
  }
}
save(CH.index,file = "./results/pbmc3k_CH.RData")