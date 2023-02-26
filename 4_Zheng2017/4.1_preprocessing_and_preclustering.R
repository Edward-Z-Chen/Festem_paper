library(Seurat)
b1_exp <- read.table("b1_exprs.txt",header = T,row.names = 1,sep = "\t")
pbmc <- CreateSeuratObject(b1_exp,project = "PBMC_2")
rm(b1_exp)
metadata <- read.table("b1_celltype.txt",header = T,sep = "\t")
pbmc <- AddMetaData(pbmc,metadata$CellType, col.name = "celltype")
rm(metadata)
saveRDS(pbmc,file = "./results/b1.rds")

# Preclustering -----------------------------------------------------------
pbmc <- readRDS("./results/b1.rds")
pbmc <- NormalizeData(pbmc)
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]
pbmc <- FindVariableFeatures(pbmc,nfeatures = 8000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc,features = VariableFeatures(pbmc),verbose = F)
# ElbowPlot(pbmc,ndims = 50)
pbmc <- FindNeighbors(pbmc,dims = 1:30)
pbmc <- FindClusters(pbmc,resolution = 0.2)
pbmc <- RunUMAP(pbmc,dims = 1:30)
DimPlot(pbmc,label = T)
DimPlot(pbmc,label = T,group.by = "celltype")
cluster.labels <- pbmc@active.ident
levels(cluster.labels) <- 1:nlevels(cluster.labels)
save(cluster.labels,file = "./results/b1_preclustering.RData")

pbmc <- FindClusters(pbmc,resolution = 0.4)
cluster.labels <- pbmc@active.ident
levels(cluster.labels) <- 1:nlevels(cluster.labels)
save(cluster.labels,file = "./results/b1_preclustering_13g.RData")

pbmc <- FindClusters(pbmc,resolution = 0.1)
cluster.labels <- pbmc@active.ident
levels(cluster.labels) <- 1:nlevels(cluster.labels)
save(cluster.labels,file = "./results/b1_preclustering_7g.RData")

# Count matrix for TN test (python) ---------------------------------------
counts <- as.matrix(pbmc@assays$RNA@counts)
save(counts,file = "./results/b1_counts_forpython.RData")

b2_exp <- read.table("b2_exprs.txt",header = T,row.names = 1,sep = "\t")
pbmc <- CreateSeuratObject(b2_exp,project = "PBMC_2")
rm(b2_exp)
metadata <- read.table("b2_celltype.txt",header = T,sep = "\t")
pbmc <- AddMetaData(pbmc,metadata$CellType, col.name = "celltype")
rm(metadata)
saveRDS(pbmc,file = "./results/b2.rds")

# Preclustering -----------------------------------------------------------
pbmc <- readRDS("./results/b2.rds")
pbmc <- NormalizeData(pbmc)
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]
pbmc <- FindVariableFeatures(pbmc,nfeatures = 8000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc,features = VariableFeatures(pbmc),verbose = F)
pbmc <- FindNeighbors(pbmc,dims = 1:30)
pbmc <- FindClusters(pbmc,resolution = 0.2)
pbmc <- RunUMAP(pbmc,dims = 1:30)
DimPlot(pbmc,label = T)
DimPlot(pbmc,label = T,group.by = "celltype")
cluster.labels <- pbmc@active.ident
levels(cluster.labels) <- 1:nlevels(cluster.labels)
save(cluster.labels,file = "./results/b2_preclustering.RData")

pbmc <- FindClusters(pbmc,resolution = 0.3)
cluster.labels <- pbmc@active.ident
levels(cluster.labels) <- 1:nlevels(cluster.labels)
save(cluster.labels,file = "./results/b2_preclustering_12g.RData")

pbmc <- FindClusters(pbmc,resolution = 0.1)
cluster.labels <- pbmc@active.ident
levels(cluster.labels) <- 1:nlevels(cluster.labels)
save(cluster.labels,file = "./results/b2_preclustering_8g.RData")

# Count matrix for TN test (python) ---------------------------------------
counts <- as.matrix(pbmc@assays$RNA@counts)
save(counts,file = "./results/b2_counts_forpython.RData")