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

# Figure 5C ---------------------------------------------------------------
library(dplyr)
library(survival)
library(survminer)
draw_survival_plot <- function(expression, metadata, title){
  genedata <- metadata[,c("Sample","Vital","Survival")]
  genedata <- cbind(genedata,FAM3C_DUSP4 = rowSums(t(expression[c("FAM3C","DUSP4"),metadata$Sample])))
  genedata <- cbind(genedata,FAM3C_DUSP4_level = genedata$FAM3C_DUSP4>=mean(genedata$FAM3C_DUSP4))
  
  genedata <- cbind(genedata,CX3CR1 = t(expression["CX3CR1",metadata$Sample]))
  genedata <- cbind(genedata,CX3CR1_level = genedata$CX3CR1>=mean(genedata$CX3CR1))
  genedata <- cbind(genedata,group = NA)
  for (i in 1:nrow(genedata)){
    if (genedata$FAM3C_DUSP4_level[i]==T & genedata$CX3CR1_level[i]==T){
      genedata$group[i] <- 1
    } else if (genedata$FAM3C_DUSP4_level[i]==T & genedata$CX3CR1_level[i]==F){
      genedata$group[i] <- 2
    } else if (genedata$FAM3C_DUSP4_level[i]==F & genedata$CX3CR1_level[i]==T){
      genedata$group[i] <- 3
    } else {
      genedata$group[i] <- 4
    }
  }
  genedata$group <- factor(genedata$group)
  result <- survival::survfit(Surv(Survival, Vital)~group, data=genedata)
  p <- survminer::ggsurvplot(result,data = genedata, conf.int=FALSE, pval=TRUE, risk.table=F, 
                   legend.labs=c("emrahi exhi","emralow exhi","emrahi exlow", "emralow exlow"), legend.title="",  
                   palette=c("#66BFBF","#E7B800","orchid2", "dodgerblue2"),
                   title=title,legend = "bottom")
  return(p)
}

expression <- read.table("./iCCA_survival/cancerdiscov_expr_matrix.txt",header = T)
metadata <- read.table("./iCCA_survival/cancerdiscov_clinical.txt",header = T,na.strings = "N/A",sep = "\t")
metadata$Sample <- paste0("X",metadata$Sample)
expression <- expression[,colnames(expression)%in%metadata$Sample]

p1 <- draw_survival_plot(expression, metadata, "Jusakul et al. (2017)")

expression <- read.table("./iCCA_survival/cptac_expr_matrix.txt",header = T)
metadata <- read.table("./iCCA_survival/cptac_clinical.txt",header = T,na.strings = "NA",sep = "\t")
metadata$Sample <- paste0("X",metadata$Sample)
metadata <- metadata[metadata$Sample%in%colnames(expression),]

p2 <- draw_survival_plot(expression, metadata, "Dong et al. (2017)")

pdf("./figures/Figure5C.pdf",width = 13.3,height = 6.6)
arrange_ggsurvplots(list(p1,p2),nrow = 1,ncol = 2)
dev.off()

# Figure 5E ---------------------------------------------------------------
expression <- read.table("./iCCA_survival/cptac_expr_matrix.txt",header = T)
metadata <- read.table("./iCCA_survival/cptac_clinical.txt",header = T,na.strings = "NA",sep = "\t")
metadata$Sample <- paste0("X",metadata$Sample)
metadata <- metadata[metadata$Sample%in%colnames(expression),]

genedata <- metadata[,c("Sample","Vital","Survival")]
genedata <- cbind(genedata,FAM3C_DUSP4 = rowSums(t(expression[c("FAM3C","DUSP4"),metadata$Sample])))
genedata <- cbind(genedata,FAM3C_DUSP4_level = genedata$FAM3C_DUSP4>=mean(genedata$FAM3C_DUSP4))

genedata <- cbind(genedata,CX3CR1 = t(expression["CX3CR1",metadata$Sample]))
genedata <- cbind(genedata,CX3CR1_level = genedata$CX3CR1>=mean(genedata$CX3CR1))
genedata <- cbind(genedata,group = NA)
for (i in 1:nrow(genedata)){
  if (genedata$FAM3C_DUSP4_level[i]==T & genedata$CX3CR1_level[i]==T){
    genedata$group[i] <- 1
  } else if (genedata$FAM3C_DUSP4_level[i]==T & genedata$CX3CR1_level[i]==F){
    genedata$group[i] <- 2
  } else if (genedata$FAM3C_DUSP4_level[i]==F & genedata$CX3CR1_level[i]==T){
    genedata$group[i] <- 3
  } else {
    genedata$group[i] <- 4
  }
}
genedata$group <- factor(genedata$group)
genedata <- cbind.data.frame(genedata,POSTN = t(expression["POSTN",]))
genedata <- cbind(genedata,POSTN_level = genedata$POSTN>=mean(genedata$POSTN))

result <- survfit(Surv(Survival, Vital)~POSTN_level, data=genedata[genedata$group==3,])
p1 <- ggsurvplot(result, conf.int=F, pval=TRUE, risk.table=F, 
                 legend.labs=c("Low","High"), legend.title="POSTN",  
                 palette=c("orchid2","#E7B800"),
                 title="emrahiexlow in Dong et al. (2017)",legend = "bottom")

result <- survfit(Surv(Survival, Vital)~POSTN_level, data=genedata[genedata$group==1,])
p2 <- ggsurvplot(result, conf.int=F, pval=TRUE, risk.table=F, 
                 legend.labs=c("Low","High"), legend.title="POSTN",  
                 palette=c("orchid2","#E7B800"),
                 title="emrahiexhi in Dong et al. (2017)",legend = "bottom")

pdf("./figures/Figure5E.pdf",width = 13.3,height = 6.6)
arrange_ggsurvplots(list(p1,p2),nrow = 1,ncol = 2)
dev.off()


# Figure S14 --------------------------------------------------------------
nonepi <- readRDS("NonEpi.rds")
nonepi <- NormalizeData(nonepi)
marker_list <- c("GNLY", "NKG7", "CCL5", "CD247","GZMB", # CD8 T
                 "EOMES","GZMK","GZMA","CCL5","GZMH","IL32","CCL4","HLA-DRB1","HLA-DPA1","COTL1","HLA-DQA1", # GZMK Tem
                 #"CMC1","TRAT1","DTHD1","PIK3R1","CRTAM", # GZMK early Tem
                 "SLC4A10","KLRB1","CCR6","DPP4","LTK", "NCR3",#MAIT https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8363247/  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6506535/
                 "IL7R","CCR7","CD27","SELL","TCF7","CD3E", # CD4 T
                 "NKG7", "GNLY", "CD247","GZMB", #NK
                 "NCAM1", # CD56dim NK (CD56=NCAM1)
                 "TNFSF10","KLRC1", # CD56bright NK
                 "CX3CR1","KLF2","ZEB2","FGFBP2","FCGR3A","GNLY","C1orf21", # Temra
                 "MS4A1", "CD79A", "CD37", # B
                 "CD68","CD14","CD163","ITGAM","ITGAX","IFIT1", # Macrophage
                 "CD14", "LYZ", "S100A9", #CD14 monocyte
                 "KLRB1","KRT81","TNFRSF18","KRT86","KLRC2","AREG", # Innate lymphoid cells https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8294935/ (ILC)
                 "IL2RA","CD4","FOXP3",# T reg
                 "CLDN5", "PECAM1", "CD34", "FLT1", "VWF", "ENG", "CDH5", # endothelial
                 "JAM2","HLA-DRB1","MCAM","SPARCL1", # blood endothelial cell
                 "CD1C","ITGAX","ITGAM","SIRPA","LILRA2","CLEC4A","CLEC10A", # myeloid cDC2
                 "EGR1","CDC14A","FOSB", #FOSB+ T
                 "SDC1","CD38","CD27","MZB1","SDC1","CD79A","IGHG1","IGKC","JCHAIN", # plasma
                 "CCL3","CCL20","STMN1","FABP5","TNFRSF9","TNFRSF18","FAM3C","DUSP4", # Terminal Tex
                 "STAT1","IRF7","IFI44L","MX1","IFI6","ISG15", # ISG+ CD8 T
                 "POSTN","COL1A2", "FAP", "DCN", "COL3A1", "COL6A1", "COL1A1", # Fibroblast
                 "CLEC9A","CADM1","XCR1","BTLA","DPP4", # Myeloid cDC1
                 "MKI67","TOP2A", # Cycling cell
                 "FCGR3A","CD14","CD86", "S100A10", "IFI30","HLA-DRB5","GBP4", # intermediate monocyte https://www.frontiersin.org/articles/10.3389/fimmu.2019.02035/full (reference 3-5)
                 "PROX1","COLEC12","CCL21" # lymphatic endothelial cell
)
marker_list <- unique(marker_list)
label.tmp <- label.list[[1]]
levels(label.tmp) <- rank(-summary(label.tmp),ties.method = "first")
levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),1:22,c("GZMK Tem","MAIT","CD4 T",
                                                              "CD56dim NK","CD56bright NK","Temra",
                                                              "B","Macrophage","Classical Monocyte",
                                                              "Innate lymphoid cells",
                                                              "Treg","Blood Endothelial Cell",
                                                              "Myeloid cDC2","FOSB+ T","Plasma",
                                                              "Terminal Tex",
                                                              "ISG+ CD8 T-like",
                                                              "Fibroblast","Myeloid cDC1","Cycling Cell",
                                                              "Intermediate Monocyte","Lymphatic Endothelial cell"))
label.tmp <- factor(label.tmp,levels = rev(levels(label.tmp)))
nonepi@active.ident <- label.tmp
nonepi <- ScaleData(nonepi,features = marker_list)

pdf("./figures/FigureS14.pdf",width = 24,height = 8)
DotPlot(nonepi, features = marker_list, cols = c("blue", "red"), 
        dot.scale = 8,idents = levels(label.tmp))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

# Figure S15 --------------------------------------------------------------
load("./results/iCCA_clustering_UMAP.RData")
for (i in 1:6){
  levels(label.list[[i]]) <- 1:nlevels(label.list[[i]])
  my.color <- hue_pal()(nlevels(label.list[[i]]))
  names(my.color) <- (1:nlevels(label.list[[i]]))
  umap.tmp <- umap.for.plot(umap.list[[1]],label.list[[i]])
  
  class_avg <- umap.tmp %>%
    group_by(cluster) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  
  plots.list[[i]] <- ggplot(umap.tmp, aes(x=UMAP_1, y=UMAP_2, color=cluster)) + 
    geom_scattermore(pointsize = 0.5) + theme_pubr()+
    theme(legend.position="none",text = element_text(size = 10)) +
    geom_text(aes(x=UMAP_1,y = UMAP_2,label = cluster), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 3)+
    labs(title = names(label.list)[i])+scale_color_manual(values = my.color,name = "Clusters")
}

pdf("./figures/FigureS15.pdf",width = 6,height = 9)
ggarrange(plotlist = plots.list,ncol = 2,nrow = 3,align = "hv")
dev.off()