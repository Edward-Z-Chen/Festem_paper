library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(scales)
library(scattermore)
if (!file.exists("./figures")){
  dir.create("./figures")
}

my.color <- hue_pal()(22)
names(my.color) <- 1:22
my.color[c(3,9)] <- my.color[c(9,3)]
umap.for.plot <- function(umap,cluster){
  umap <- umap@cell.embeddings
  umap <- as.data.frame(umap)
  cbind(umap,cluster = cluster)
}

# Figure 5A ---------------------------------------------------------------
load("./results/iCCA_clustering_UMAP.RData")
label.tmp <- label.list[[1]]
levels(label.tmp) <- rank(-summary(label.tmp),ties.method = "first")
umap.tmp <- umap.for.plot(umap.list[[1]],label.tmp)

class_avg <- umap.tmp %>%
  group_by(cluster) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
class_avg <- cbind(class_avg,class_avg[,1])
colnames(class_avg)[4] <- "cluster.anno"
class_avg[,4] <- plyr::mapvalues(class_avg[,4],1:22,c("GZMK Tem","MAIT","CD4 T",
                                                      "CD56Dim NK","CD56Bright NK","Temra",
                                                      "B","Macrophage","Classical Monocyte",
                                                      "Innate lymphoid cells",
                                                      "Treg","Blood Endothelial Cell",
                                                      "Myeloid cDC2","FOSB+ T","Plasma",
                                                      "Terminal Tex",
                                                      "ISG+ CD8 T-like",
                                                      "Fibroblast","Myeloid cDC1","Cycling Cell",
                                                      "Intermediate Monocyte","Lymphatic Endothelial cell"))
pdf("./figures/Figure5A.pdf",width = 4, height = 4)
ggplot(umap.tmp, aes(x=UMAP_1, y=UMAP_2, color=cluster)) + 
  geom_scattermore(pointsize = 0.5) + theme_pubr()+
  theme(legend.position="none",text = element_text(size = 15)) +
  # ggrepel::geom_label_repel(data = class_avg,aes(x=UMAP_1,y = UMAP_2,label = cluster.anno),label.padding = unit(0.1,'cm'),label.size = 0.15)+
  geom_text(aes(x=UMAP_1,y = UMAP_2,label = cluster.anno), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 3)+
  labs(title = "Festem")+scale_color_manual(values = my.color,name = "Clusters")
dev.off()


## zoom in ------------------------------------------------------------
umap.tmp <- filter(umap.tmp,cluster %in% c(1:3,6,10,11,14,16,17))
umap.tmp <- filter(umap.tmp,UMAP_1>-5 | UMAP_2<0)
umap.tmp <- filter(umap.tmp,UMAP_2<8)
umap.tmp <- filter(umap.tmp,UMAP_1>-5)
class_avg <- umap.tmp %>%
  group_by(cluster) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
class_avg <- cbind(class_avg,class_avg[,1])
colnames(class_avg)[4] <- "cluster.anno"
class_avg[,4] <- plyr::mapvalues(class_avg[,4],1:22,c("GZMK Tem","MAIT","CD4 T",
                                                      "CD56Dim NK","CD56Bright NK","Temra",
                                                      "B","Macrophage","Classical Monocyte",
                                                      "Innate lymphoid cells",
                                                      "Treg","Blood Endothelial Cell",
                                                      "Myeloid cDC2","FOSB+ T","Plasma",
                                                      "Terminal Tex",
                                                      "ISG+ CD8 T-like",
                                                      "Fibroblast","Myeloid cDC1","Cycling Cell",
                                                      "Intermediate Monocyte","Lymphatic Endothelial cell"))

pdf("./figures/Figure5A_zoomin.pdf",width = 4, height = 4)
ggplot(umap.tmp, 
       aes(x=UMAP_1, y=UMAP_2, color=cluster)) + 
  geom_scattermore(pointsize = 0.5) + theme_pubr()+
  theme(legend.position="none",text = element_text(size = 15)) +
  # ggrepel::geom_label_repel(data = class_avg,aes(x=UMAP_1,y = UMAP_2,label = cluster.anno),label.padding = unit(0.1,'cm'),label.size = 0.15)+
  geom_text(aes(x=UMAP_1,y = UMAP_2,label = cluster.anno), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 3)+
  labs(title = "Festem")+scale_color_manual(values = my.color,name = "Clusters")
dev.off()

# Figure 5B ---------------------------------------------------------------
library(cowplot)
library(patchwork)
load("./results/iCCA_clustering_UMAP.RData")
gene.names <- gene.list[[1]]
nonepi <- readRDS("NonEpi.rds")
nonepi <- NormalizeData(nonepi)
nonepi <- nonepi[,rownames(metadata)]
nonepi@active.ident <- factor(label.list[[1]])
names(nonepi@active.ident) <- colnames(nonepi)
marker_list <- c("CX3CR1","ZEB2","GZMH","FGFBP2",
                 "FCGR3A","PRSS23","GNLY","NKG7","KLRD1",
                 "PLEK","C1orf21",
                 "GZMB",
                 "CCL3","CCL20","STMN1","FABP5",
                 "TNFRSF9","TNFRSF18",
                 "FAM3C","DUSP4"
)
nonepi <- ScaleData(nonepi,marker_list)
marker_exp <- data.frame(exp = as.numeric(nonepi@assays$RNA@scale.data[marker_list,plots.list[[1]]$data$cluster=="5"]),
                         gene = rep(marker_list,sum(plots.list[[1]]$data$cluster=="5")),
                         id = "Temra")
marker_exp <- rbind(marker_exp,data.frame(exp = as.numeric(nonepi@assays$RNA@scale.data[marker_list,plots.list[[1]]$data$cluster=="0"]),
                                          gene = rep(marker_list,sum(plots.list[[1]]$data$cluster=="0")),
                                          id = "GZMK Tem"))
marker_exp <- rbind(marker_exp,data.frame(exp = as.numeric(nonepi@assays$RNA@scale.data[marker_list,plots.list[[1]]$data$cluster=="2"]),
                                          gene = rep(marker_list,sum(plots.list[[1]]$data$cluster=="2")),
                                          id = "CD4 T"))
marker_exp <- rbind(marker_exp,data.frame(exp = as.numeric(nonepi@assays$RNA@scale.data[marker_list,plots.list[[1]]$data$cluster=="10"]),
                                          gene = rep(marker_list,sum(plots.list[[1]]$data$cluster=="10")),
                                          id = "Treg"))
marker_exp <- rbind(marker_exp,data.frame(exp = as.numeric(nonepi@assays$RNA@scale.data[marker_list,plots.list[[1]]$data$cluster=="15"]),
                                          gene = rep(marker_list,sum(plots.list[[1]]$data$cluster=="15")),
                                          id = "Terminal Tex"))
set.seed(321)
noise <- rnorm(n = nrow(marker_exp)) / 100000
marker_exp$exp <- marker_exp$exp+noise
marker_exp$gene <- factor(marker_exp$gene,levels = marker_list)

my.color <- hue_pal()(22)
names(my.color) <- 1:22
my.color[c(3,9)] <- my.color[c(9,3)]
names(my.color) <- plyr::mapvalues(names(my.color),1:22,c("GZMK Tem","MAIT","CD4 T",
                                                          "CD56Dim NK","CD56Bright NK","Temra",
                                                          "B","Macrophage","Classical Monocyte",
                                                          "Innate lymphoid cells",
                                                          "Treg","Blood Endothelial Cell",
                                                          "Myeloid cDC2","FOSB+ T","Plasma",
                                                          "Terminal Tex",
                                                          "ISG+ CD8 T-like",
                                                          "Fibroblast","Myeloid cDC1","Cycling Cell",
                                                          "Intermediate Monocyte","Lymphatic Endothelial cell"))
color <- my.color[unique(marker_exp$id)]

pdf("./figures/Figure5B.pdf",width = 4, height = 6)
ggplot(marker_exp, aes(x = factor(id,levels = unique(id)),
                       y = exp, fill = id)) +
  geom_violin(scale = "width", adjust = 0.8, trim = T) +
  scale_y_continuous(expand = c(0, 0), position="left", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(gene), scales = "free") +
  cowplot::theme_cowplot(font_size = 12) +
  scale_fill_discrete(type = color)+
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.right = element_text(angle = 0),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()) +
  ggtitle("")+xlab("") + ylab("Expression Level")
dev.off()
