library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(scales)
library(RColorBrewer)
my.color <- hue_pal()(13)
names(my.color) <- 1:13
umap.for.plot <- function(umap,cluster){
  umap <- umap@cell.embeddings
  umap <- as.data.frame(umap)
  cbind(umap,cluster = cluster)
}
if (!file.exists("./figures")){
  dir.create("./figures")
}
# Figure 3 (A) -- left top and Figure S5 top ------------------------------------------------
## Summarising results from various DEG detection methods
results <- matrix(nrow = nrow(pbmc),ncol = 19,
                  dimnames = list(rownames(pbmc),
                                  c("Festem","Festem_gamma0.01","Festem_7g","Festem_10g",
                                    "DEseq2", "DEseq2-full",       
                                    "EdgeR", "EdgeR-full","MAST",              
                                    "Wilcoxon","MAST-f",            
                                    "Wilcoxon-f","FC",
                                    "singleCellHaystack-PCA",
                                    "singleCellHaystack-UMAP",
                                    "singleCellHaystack-PCA-10pc",
                                    "singleCellHaystack-UMAP-10pc","TN_test","ROSeq")))
pbmc <- readRDS("./results/pbmc3k.rds")
pbmc <- as.SingleCellExperiment(pbmc)
load("./results/pbmc3k_label.RData")
cluster.label[c("CGGGCATGACCCAA-1","CTTGATTGATCTTC-1")] <- "Platelet"
pbmc <- pbmc[,cluster.label!="Platelet"]
cluster.label <- factor(cluster.label[cluster.label!="Platelet"])
pbmc <- AddMetaData(pbmc,cluster.label,"HVG")
FC.list <- vector("list",8)
for (i in 1:8){
  FC.list[[i]] <- FoldChange(pbmc,ident.1 = levels(cluster.label)[i],
                             group.by = "HVG")[,1]
}
FC.list <- matrix(unlist(FC.list),ncol = 8,byrow = F)
FC <- apply(FC.list,1,function(x){max(abs(x))})
results[,"FC"] <- FC
load("./results/pbmc_Festem.RData")
results[colnames(em.result),"Festem"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result9)[em.result9[2,]<0],"Festem"] <- NA
load("./results/pbmc_Festem_gamma0.01.RData")
results[colnames(em.result),"Festem_gamma0.01"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result9)[em.result9[2,]<0],"Festem_gamma0.01"] <- NA
load("./results/pbmc_Festem_7g.RData")
results[colnames(em.result),"Festem_7g"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result9)[em.result9[2,]<0],"Festem_7g"] <- NA
load("./results/pbmc_Festem_10g.RData")
results[colnames(em.result),"Festem_10g"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result9)[em.result9[2,]<0],"Festem_10g"] <- NA
load("./results/pbmc_DEseq2.RData")
results[DEseq.results@rownames,"DESeq2"] <- DEseq.results@listData$padj
results[DEseq.results@rownames,"DESeq2-full"] <- p.adjust(DEseq.results@listData$pvalue,"BH")
load("./results/pbmc_edgeR.RData")
results[names(EdgeR.result),"EdgeR"] <- p.adjust(EdgeR.result,"BH")
load("./results/pbmc_edgeR_full.RData")
results[names(EdgeR.result),"EdgeR-full"] <- p.adjust(EdgeR.result,"BH")
load("./results/pbmc_mast.RData")
results[rownames(mast.results),"MAST"] <- p.adjust(mast.results[,3],"BH")
load("./results/pbmc_wilcoxon.RData")
wil.result <- matrix(unlist(wil.result),ncol = 8,byrow = F)
wil.result <- apply(wil.result,1,function(x){min(x)*8})
results[,"Wilcoxon"] <- p.adjust(wil.result,"BH")
results[,"Wilcoxon-f"] <- results[,"Wilcoxon"]
results[results[,"FC"]<=0.2,"Wilcoxon-f"] <- NA
results[,"MAST-f"] <- results[,"MAST"]
results[results[,"FC"]<=0.2,"MAST-f"] <- NA
load("./results/pbmc_haystack.RData")
results[,"singleCellHaystack-PCA"] <- exp(haystack_pca$results$log.p.adj)
results[,"singleCellHaystack-UMAP"] <- exp(haystack_umap$results$log.p.adj)
load("./results/pbmc_ROSeq.RData")
results[,"ROSeq"] <- p.adjust(roseq.tmp,"BH")

tn_test <- read.csv("./results/pbmc_TN_test.csv",
                    header = F)
tn_test <- as.matrix(tn_test)
tn_test <- apply(tn_test,1,function(x){min(x)*length(x)})
results[,"TN_test"] <- p.adjust(tn_test,"BH")
save(results,file = "./results/pbmc3k_DEG_results.RData")

## Figure S5 top left -----------------------------------------------------------
calc_house_percent <- function(DEG_list,house_list){
  sum(house_list%in%DEG_list)/length(DEG_list)
}
calc_house_retrieve <- function(DEG_list,house_list){
  sum(DEG_list%in%house_list)/length(house_list)
}
load("./results/pbmc3k_silver_standard.RData")
g1 <- rownames(moran_h)[abs(moran_h[,1])<0.02 & moran_h[,2]>0.05]
moran_nonh <- moran_nonh[!is.na(moran_nonh[,1]),]
g2 <- rownames(moran_nonh)[moran_nonh[,1]>0.1]
g2 <- g2[gene_percent[g2]>0.05]
p_values <- results[,-13]
p_values_full <- p_values
p_values <- p_values[c(g1,g2),]

t_list <- seq(0.01,0.1,0.01)
power_frame <- data.frame()
for (i in 1:length(t_list)){
  for (j in 1:ncol(p_values)){
    reject_gene <- rownames(p_values)[p_values[,j]<t_list[i] & (!is.na(p_values[,j]))]
    power_frame <- rbind(power_frame,
                         data.frame(method = colnames(p_values)[j],
                                    cutoff = t_list[i],
                                    precision = 1-calc_house_percent(reject_gene[reject_gene%in%c(g1,g2)],g1),
                                    recall = calc_house_retrieve(reject_gene,g2),
                                    number = sum(p_values_full[,j]<t_list[i] & (!is.na(p_values_full[,j])))))
  }
}
power_frame <- cbind(power_frame,
                     F_score = 2/(1/power_frame$precision + 1/power_frame$recall))
power_frame2 <- filter(power_frame,cutoff == 0.05)

double_exp <- function() { 
  trans <- function(x) {    
    exp(exp(x))
  }
  inv <- function(x) {
    log(log(x))
  }
  
  # return the transformation
  return(trans_new("double_exp", trans, inv))
}

my_color <- c(rep("#DE2D26",4),rep("#1F78B4",2),rep("#33A02C",2),"#FF7F00","#6A3D9A",
              "#FF7F00","#6A3D9A","#FF99D7",rep("#B15928",4),"grey50")
names(my_color) <- c("Festem","Festem_gamma0.01","Festem_7g","Festem_10g",
                     "DEseq2","DEseq2-full","EdgeR","EdgeR-full",
                     "MAST","Wilcoxon","MAST-f","Wilcoxon-f","ROSeq",
                     "singleCellHaystack-PCA","singleCellHaystack-UMAP",
                     "singleCellHaystack-PCA-10pc","singleCellHaystack-UMAP-10pc","TN_test")
my_shape <- c(16,15,17,18,16,7,16,7,16,16,9,9,16,10,16,12,13,16)
names(my_shape) <- names(my_color)

pdf("./figures/FigureS5_pbmc3k_left.pdf",height = 6,width = 6)
ggplot(data = power_frame2,mapping = aes(x = recall, y = precision))+
  geom_point(aes(color = method, shape = method),
             size = 4,stroke = 1)+
  scale_x_continuous(limits = c(0.1,1))+
  scale_y_continuous(limits = c(0.2,1))+
  geom_function(fun = ~ {1/(2/0.9-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.8-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.7-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.6-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.5-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.4-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.3-1/.x)},color = "grey75",linetype = "dashed")+
  coord_trans(y=double_exp())+
  theme_pubr()+
  geom_hline(yintercept = 1-0.05,linetype = "dashed",color = "red")+
  scale_color_manual(values = my_color)+
  scale_shape_manual(values = my_shape)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 4))
dev.off()

## Figure S5 top right -----------------------------------------------------------
num_reject <- apply(p_values_full,2,function(x){sum(x<0.05,na.rm = T)})
num_reject <- data.frame(method = colnames(p_values_full),
                         num = num_reject)
num_reject <- num_reject[order(num_reject$num,decreasing = F),]
num_reject$method <- factor(num_reject$method,levels = num_reject$method)

pdf("./figures/FigureS5_pbmc3k_right.pdf",height = 8,width = 8)
ggplot(num_reject,aes(x = method,y = num,fill = method))+
  geom_bar(stat = "identity")+
  theme_pubr()+
  geom_text(aes(label=num), vjust=1.6, color="white", size=3.5)+
  scale_fill_manual(values = my_color)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90))+
  guides(color = guide_legend(nrow = 4))
dev.off()

## Figure 3 (A) -- top left ------------------------------------------------
power_frame2 <- power_frame2[c(1,5,7,9,10,14,17,18),]
power_frame2[c(1,6),1] <- c("Festem","singleCellHaystack")
my_color <- my_color[c(1,5,7,9,10,13,15,18)]
names(my_color)[7] <- "singleCellHaystack"
my_shape <- c(16,15,17,18,7,8,9,10)
names(my_shape) <- names(my_color)

pdf("./figures/Figure3A_topleft.pdf",height = 4,width = 3.5)
ggplot(data = power_frame2,mapping = aes(x = recall, y = precision))+
  geom_point(aes(color = method, shape = method),
             size = 4,stroke = 1)+
  scale_x_continuous(limits = c(0.1,1))+
  scale_y_continuous(limits = c(0.2,1))+
  coord_trans(y=double_exp())+
  geom_function(fun = ~ {1/(2/0.9-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.8-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.7-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.6-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.5-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.4-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.3-1/.x)},color = "grey75",linetype = "dashed")+
  theme_pubr()+
  geom_hline(yintercept = 1-0.05,linetype = "dashed",color = "red")+
  scale_color_manual(values = my_color)+
  scale_shape_manual(values = my_shape)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 2))
dev.off()

num_reject <- num_reject[num_reject$method %in% c(power_frame2$method,"singleCellHaystack-UMAP"),]
num_reject$method <- factor(num_reject$method)
levels(num_reject$method) <- plyr::mapvalues(levels(num_reject$method),
                                             c("singleCellHaystack-UMAP"),
                                             c("singleCellHaystack"))

pdf("./figures/Figure3A_topleft_inset.pdf",height = 4,width = 4)
ggplot(num_reject,aes(x = method,y = num,fill = method))+
  geom_bar(stat = "identity",width = 1)+
  theme_pubr()+
  geom_text(aes(label=num), vjust=1.6, color="white", size=3.5)+
  scale_fill_manual(values = my_color)+
  theme(legend.position = "bottom",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  guides(color = guide_legend(nrow = 1))+
  labs(x = NULL,y = "Number")
dev.off()


# Figure 4 ----------------------------------------------------------------
load("./results/pbmc3k_clustering_UMAP.RData")
my.color <- hue_pal()(13)
names(my.color) <- 1:13
## Figure 4A --left ----------------------------------------------------------------
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
class_avg[,4] <- plyr::mapvalues(class_avg[,4],1:10,c("Memory CD4 T","Naive CD4 T","CD14+ Monocyte",
                                                      "B","CD8 T","FCGR3A+ Monocyte","NK","M-MDSC-like",
                                                      "CD27-CD4+\n Memory T","DC"))
pdf("./figures/Figure4A_left.pdf",width = 4,height = 4)
ggplot(umap.tmp, aes(x=UMAP_1, y=UMAP_2, color=cluster)) + 
  geom_point(cex=0.1) + theme_pubr()+
  theme(legend.position="none",text = element_text(size = 15)) +
  # ggrepel::geom_label_repel(data = class_avg,aes(x=UMAP_1,y = UMAP_2,label = cluster.anno),label.padding = unit(0.1,'cm'),label.size = 0.15)+
  geom_text(aes(x=UMAP_1,y = UMAP_2,label = cluster.anno), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 3)+
  labs(title = names(label.list)[1])+scale_color_manual(values = my.color,name = "Clusters")
dev.off()

## Competing methods ----------------------------------------------------------------
for (i in 1:length(label.list)){
  label.tmp <- label.list[[i]]
  levels(label.tmp) <- rank(-summary(label.tmp),ties.method = "first")
  if (i == 2) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),3:9,c(4,5,8,6,3,7,10))
  }
  if (i == 3) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),2:8,c(3,2,4,5,6,7,10))
  }
  if (i == 4) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),5:9,c(9,7,6,5,10))
  }
  if (i == 5) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),3:9,c(4,8,5,3,9,7,6))
  }
  if (i == 6) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),2:12,c(3,4,2,6,9,7,5,11,12,10,13))
  }
  umap.tmp <- umap.for.plot(umap.list[[i]],label.tmp)
  
  class_avg <- umap.tmp %>%
    group_by(cluster) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  plots.list[[i]] <- ggplot(umap.tmp, aes(x=UMAP_1, y=UMAP_2, color=cluster)) + 
    geom_point(cex=0.1) + theme_pubr()+theme(legend.position="none",text = element_text(size = 15)) +
    # ggrepel::geom_label_repel(data = class_avg,aes(x=UMAP_1,y = UMAP_2,label = allgenes),label.padding = unit(0.1,'cm'))
    geom_text(aes(x=UMAP_1,y = UMAP_2,label = cluster), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 4)+
    labs(title = names(label.list)[i])+scale_color_manual(values = my.color,name = "Clusters")
}
pdf("./figures/umaps.pdf",height = 6,width = 10)
ggarrange(plotlist = plots.list,nrow = 2,ncol = 3,legend = "none")
dev.off()

## Figure 4A --right ----------------------------------------------------------------
for (i in 1:length(label.list)){
  label.tmp <- label.list[[i]]
  levels(label.tmp) <- rank(-summary(label.tmp),ties.method = "first")
  if (i == 2) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),3:9,c(4,5,8,6,3,7,10))
  }
  if (i == 3) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),2:8,c(3,2,4,5,6,7,10))
  }
  if (i == 4) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),5:9,c(9,7,6,5,10))
  }
  if (i == 5) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),3:9,c(4,8,5,3,9,7,6))
  }
  if (i == 6) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),2:12,c(3,4,2,6,9,7,5,11,12,10,13))
  }
  umap.tmp <- umap.for.plot(umap.list[[1]],label.tmp)
  
  class_avg <- umap.tmp %>%
    group_by(cluster) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  umap.tmp <- filter(umap.tmp, UMAP_1<0 & UMAP_2>-5)
  umap.tmp$cluster <- factor(umap.tmp$cluster)
  class_avg <- class_avg[class_avg$cluster %in% levels(umap.tmp$cluster),]
  plots.list[[i]] <- ggplot(umap.tmp, aes(x=UMAP_1, y=UMAP_2, color=cluster)) + 
    geom_point(cex=0.1) + theme_pubr()+theme(legend.position="none",text = element_text(size = 15)) +
    # ggrepel::geom_label_repel(data = class_avg,aes(x=UMAP_1,y = UMAP_2,label = allgenes),label.padding = unit(0.1,'cm'))
    geom_text(aes(x=UMAP_1,y = UMAP_2,label = cluster), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 4)+
    labs(title = names(label.list)[i])+
    scale_x_continuous(limits = c(-10,-3))+
    scale_y_continuous(limits = c(-5,8))
  if (i %in% c(2,3)){
    plots.list[[i]] <- plots.list[[i]]+scale_color_manual(values = my.color[c(5)],name = "Clusters")
  } else{
    plots.list[[i]] <- plots.list[[i]]+scale_color_manual(values = my.color[c(9)],name = "Clusters")
  }
}
pdf("./figures/Figure4A_right.pdf",width = 5,height = 4)
ggarrange(plotlist = plots.list,
          nrow = 2,ncol = 3,legend = "none")
dev.off()


# Figure 4B ---------------------------------------------------------------
library(cowplot)
pbmc <- readRDS("./results/pbmc3k.rds")
ref <- pbmc@active.ident
ref[c("CGGGCATGACCCAA-1","CTTGATTGATCTTC-1")] <- "Platelet"
pbmc <- pbmc[,ref!="Platelet"]
load("./results/pbmc3k_clustering_UMAP.RData")
marker_list <- c("CD27","SELL","CCR7","MAL","LEF1","CCL5",
                 "CST7","NKG7","S100A11","IL7R","S100A4")
marker_exp <- data.frame(exp = as.numeric(pbmc@assays$RNA@data[marker_list,plots.list[[1]]$data$cluster=="8"]),
                         gene = rep(marker_list,sum(plots.list[[1]]$data$cluster=="8")),
                         id = "Festem")
marker_exp <- rbind(marker_exp,data.frame(exp = as.numeric(pbmc@assays$RNA@data[marker_list,plots.list[[5]]$data$cluster=="6"]),
                                          gene = rep(marker_list,sum(plots.list[[5]]$data$cluster=="6")),
                                          id = "devianceFS"))
marker_exp <- rbind(marker_exp,data.frame(exp = as.numeric(pbmc@assays$RNA@data[marker_list,plots.list[[6]]$data$cluster=="5"]),
                                          gene = rep(marker_list,sum(plots.list[[6]]$data$cluster=="5")),
                                          id = "trendVar"))
marker_exp <- rbind(marker_exp,data.frame(exp = as.numeric(pbmc@assays$RNA@data[marker_list,plots.list[[4]]$data$cluster=="4"]),
                                          gene = rep(marker_list,sum(plots.list[[4]]$data$cluster=="4")),
                                          id = "DUBStepR"))
set.seed(321)
noise <- rnorm(n = nrow(marker_exp)) / 100000
marker_exp$exp <- marker_exp$exp+noise
marker_exp$gene <- factor(marker_exp$gene,levels = marker_list)
my.color <- c("#DE2D26","#4DAF4A","#FF7F00","#C0B236")
names(my.color) <- c("Festem","DUBStepR","devianceFS","trendVar")

pdf("./figures/Figure4B.pdf",height = 5.5,width = 4)
ggplot(marker_exp, aes(x = factor(id,levels = c("Festem","devianceFS","trendVar","DUBStepR")),
                       y = exp, fill = id)) +
  geom_violin(scale = "width", adjust = 1, trim = T) +
  scale_y_continuous(expand = c(0, 0), position="left", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(gene), scales = "free") +
  theme_cowplot(font_size = 12) +
  scale_fill_discrete(type = my.color)+
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.right = element_text(angle = 0),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()) +
  ggtitle("CD27-CD4+ Memory T (Cluster 9)")+xlab("") + ylab("Expression Level")
dev.off()


# Figure 4C & D ---------------------------------------------------------------
library(ggalluvial)
library(tidyverse)
library(scales)
my.color <- hue_pal()(13)
names(my.color) <- 1:13
names(my.color) <- plyr::mapvalues(names(my.color),1:10,c("Memory CD4 T","Naive CD4 T","CD14+ Monocyte",
                                                          "B","CD8 T","FCGR3A+ Monocyte","NK","M-MDSC-like",
                                                          "CD27-CD4+\n Memory T","DC"))
load("./results/pbmc3k_clustering_UMAP.RData")
for (i in 1:length(label.list)){
  label.tmp <- label.list[[i]]
  levels(label.tmp) <- rank(-summary(label.tmp),ties.method = "first")
  if (i == 2) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),3:9,c(4,5,8,6,3,7,10))
  }
  if (i == 3) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),2:8,c(3,2,4,5,6,7,10))
  }
  if (i == 4) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),5:9,c(9,7,6,5,10))
  }
  if (i == 5) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),3:9,c(4,8,5,3,9,7,6))
  }
  if (i == 6) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),2:12,c(3,4,2,6,9,7,5,11,12,10,13))
  }
  umap.tmp <- umap.for.plot(umap.list[[i]],label.tmp)
  
  class_avg <- umap.tmp %>%
    group_by(cluster) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  plots.list[[i]] <- ggplot(umap.tmp, aes(x=UMAP_1, y=UMAP_2, color=cluster)) + 
    geom_point(cex=0.1) + theme_pubr()+theme(legend.position="none",text = element_text(size = 15)) +
    # ggrepel::geom_label_repel(data = class_avg,aes(x=UMAP_1,y = UMAP_2,label = allgenes),label.padding = unit(0.1,'cm'))
    geom_text(aes(x=UMAP_1,y = UMAP_2,label = cluster), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 4)+
    labs(title = names(label.list)[i])+scale_color_manual(values = my.color,name = "Clusters")
}

label_data <- rownames(plots.list[[1]]$data)
label_data <- data.frame(label_data)
for (i in 1:length(label.list)){
  label_data <- cbind.data.frame(label_data,plots.list[[i]]$data$cluster)
}
label_data <- cbind.data.frame(label_data,1)
colnames(label_data) <- c("cell",names(label.list),"freq")
for (i in 2:(ncol(label_data)-1)){
  levels(label_data[,i]) <- plyr::mapvalues(levels(label_data[,i]),1:10,c("Memory CD4 T","Naive CD4 T","CD14+ Monocyte",
                                                                          "B","CD8 T","FCGR3A+ Monocyte","NK","M-MDSC-like",
                                                                          "CD27-CD4+\n Memory T","DC"))
}

pdf("./figures/Figure4C.pdf",height = 5,width = 4)
ggplot(data = label_data,
       aes(axis2 = devianceFS, axis1 = Festem, 
           y = freq)) +
  geom_alluvium(aes(fill = Festem)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Festem", "devianceFS"),
                   expand = c(0.15, 0.05)) +
  # ggfittext::geom_fit_text(stat = "stratum") +
  theme_void()+theme(legend.position = "none")+
  scale_fill_manual(values = my.color)
dev.off()

pdf("./figures/Figure4D.pdf",height = 5,width = 4)
ggplot(data = label_data,
       aes(axis2 = HVGvst, axis1 = Festem, 
           y = freq)) +
  geom_alluvium(aes(fill = Festem)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Festem", "HVGvst"),
                   expand = c(0.15, 0.05)) +
  # ggfittext::geom_fit_text(stat = "stratum") +
  theme_void()+theme(legend.position = "none")+
  scale_fill_manual(values = my.color)
dev.off()


# Figure S1 ---------------------------------------------------------------
pbmc <- readRDS("./results/pbmc3k.rds")
ref <- pbmc@active.ident
ref[c("CGGGCATGACCCAA-1","CTTGATTGATCTTC-1")] <- "Platelet"
pbmc <- pbmc[,ref!="Platelet"]
load("./results/pbmc3k_hvggenes.RData")
load("./results/pbmc3k_DEG_results.RData")
EM <- rownames(results)[results[,1]<=0.05]
EM <- EM[!is.na(EM)]
tmp <- hvgvst %in% EM
load("Housekeeping_GenesHuman.RData")
tmp <- !tmp & (hvgvst %in% Housekeeping_Genes$Gene.name)

## Figure S1 (A) ---------------------------------------------------------------
load("./results/pbmc3k_clustering_UMAP.RData")
marker_list <- hvgvst[which(tmp)[1:10]]
feature_plot_list <- vector("list",length = length(marker_list))

for (i in 1:length(marker_list)){
  data.tmp <- plots.list[[1]]$data
  data.tmp$cluster <- factor(data.tmp$cluster)
  exp.tmp <- pbmc@assays$RNA@data[marker_list[i],]
  data.tmp <- cbind.data.frame(data.tmp,exp = exp.tmp)
  feature_plot_list[[i]] <- 
    ggplot(mapping = aes(x=UMAP_1, y=UMAP_2)) + 
    geom_point(data = data.tmp,
               aes(color=exp),
               cex=0.5,shape = 16)+ 
    scale_color_gradient(low = "grey75",high = "blue")+
    theme_pubr()+theme(legend.position="none") +
    labs(title = marker_list[i])
}
pdf("./figures/FigureS1A.pdf",width = 16,height = 8)
ggarrange(plotlist = feature_plot_list,ncol = 5,nrow = 2,legend = "none")
dev.off()

## Figure S1 (B) ---------------------------------------------------------------
EM <- rownames(results)[results[,1]<=0.05]
EM <- EM[!is.na(EM)]
tmp <- EM %in% hvgvst
tmp <- !tmp & (!EM %in% Housekeeping_Genes$Gene.name)

marker_list <- c("GZMM","BLVRB","CSF3R","CD8A","NUP214","TSPO",
                 "IL2RB","ITGB2","RNF130","CSF1R")
feature_plot_list <- vector("list",length = length(marker_list))
for (i in 1:length(marker_list)){
  data.tmp <- plots.list[[1]]$data
  data.tmp$cluster <- factor(data.tmp$cluster)
  exp.tmp <- pbmc@assays$RNA@data[marker_list[i],]
  data.tmp <- cbind.data.frame(data.tmp,exp = exp.tmp)
  feature_plot_list[[i]] <- 
    ggplot(mapping = aes(x=UMAP_1, y=UMAP_2)) + 
    geom_point(data = data.tmp,
               aes(color=exp),
               cex=0.5,shape = 16)+ 
    scale_color_gradient(low = "grey75",high = "blue")+
    theme_pubr()+theme(legend.position="none") +
    labs(title = marker_list[i])
}
pdf("./figures/FigureS1B.pdf",width = 16,height = 8)
ggarrange(plotlist = feature_plot_list,ncol = 5,nrow = 2,legend = "none")
dev.off()

# Figure S4A --------------------------------------------------------------
marker_list <- c("AURKAIP1","CYC1","NAP1L4","SF3B5","RPN1","MRPL54")
feature_plot_list <- vector("list",length = length(marker_list))

for (i in 1:length(marker_list)){
  data.tmp <- plots.list[[1]]$data
  data.tmp$cluster <- factor(data.tmp$cluster)
  exp.tmp <- pbmc@assays$RNA@data[marker_list[i],]
  data.tmp <- cbind.data.frame(data.tmp,exp = exp.tmp)
  feature_plot_list[[i]] <- 
    ggplot(mapping = aes(x=UMAP_1, y=UMAP_2)) + 
    geom_point(data = data.tmp,
               aes(color=exp),
               cex=0.5,shape = 16)+ 
    scale_color_gradient(low = "grey75",high = "blue")+
    theme_pubr()+theme(legend.position="none") +
    labs(title = marker_list[i])
}
pdf("./figures/FigureS4A.pdf",width = 8,height = 6)
ggarrange(plotlist = feature_plot_list,ncol = 3,nrow = 2,legend = "none")
dev.off()

# Figure S6 -- top left ---------------------------------------------------
load("./results/pbmc3k_hvggenes.RData")
pbmc <- readRDS("./results/pbmc3k.rds")
pbmc <- FindVariableFeatures(pbmc,nfeatures = nrow(pbmc))
hvgvst <- VariableFeatures(pbmc)
pbmc <- FindVariableFeatures(pbmc,nfeatures = nrow(pbmc),selection.method = "disp")
hvgdisp <- VariableFeatures(pbmc)

load("./results/pbmc3k_silver_standard.RData")
g1 <- rownames(moran_h)[abs(moran_h[,1])<0.02 & moran_h[,2]>0.05]
moran_nonh <- moran_nonh[!is.na(moran_nonh[,1]),]
g2 <- rownames(moran_nonh)[moran_nonh[,1]>0.1]
g2 <- g2[gene_percent[g2]>0.05]
load("./results/pbmc3k_clustering_UMAP.RData")
## ROC curve (codes modified from https://gist.github.com/charly06/91578196fc615c5a79c7174318be4349#file-ggrocs-r)
#' Functions plots multiple 'roc' objects into one plot
#' @param rocs
#'   A list of 'roc' objects. Every list item has a name.
#' @param breaks
#'   A vector of integers representing ticks on the x- and y-axis
#' @param legentTitel
#'   A string which is used as legend titel
ggrocs <- function(rocs, breaks = seq(0,1,0.1), legendTitel = "Legend") {
  if (length(rocs) == 0) {
    stop("No ROC objects available in param rocs.")
  } else {
    require(plyr)
    # Store all sensitivities and specifivities in a data frame
    # which an be used in ggplot
    RocVals <- plyr::ldply(names(rocs), function(rocName) {
      if(class(rocs[[rocName]]) != "roc") {
        stop("Please provide roc object from pROC package")
      }
      data.frame(
        fpr = rev(1-rocs[[rocName]]$specificities),
        tpr = rev(rocs[[rocName]]$sensitivities),
        names = rep(rocName, length(rocs[[rocName]]$sensitivities)),
        stringAsFactors = T
      )
    })
    
    my.color <- c("#DE2D26","#6A3D9A","#9BA3EB","#4DAF4A","#FF7F00","#C0B236")
    names(my.color) <- c("Festem","HVGvst","HVGdisp","DUBStepR","devianceFS","trendVar")
    # aucAvg <- mean(sapply(rocs, "[[", "auc"))
    auc_tmp <- sapply(rocs, "[[", "auc")
    
    rocPlot <- ggplot(RocVals, aes(x = fpr, y = tpr, colour = names)) +
      # geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0), alpha = 0.5, colour = "gray") + 
      geom_step() +
      scale_x_continuous(name = "False Positive Rate (1 - Specificity)",limits = c(0,1), breaks = breaks) + 
      scale_y_continuous(name = "True Positive Rate (Sensitivity)", limits = c(0,1), breaks = breaks) +
      theme_pubr() + 
      coord_equal() + 
      #annotate("text", x = 0.1, y = 0.1, vjust = 0, label = paste("AUC =",sprintf("%.3f",aucAvg))) +
      guides(colour = guide_legend(legendTitel)) +
      theme(axis.ticks = element_line(color = "black"),legend.position = "right")+
      scale_color_manual(values = my.color,
                         labels = paste0(names(my.color)," (AUC = ",round(auc_tmp,3),")"))
    
    rocPlot
  }
}

library(pROC)
pROC_list <- vector("list",6)
set.seed(321)
gene.list <- list(Festem = c(EM,setdiff(rownames(pbmc),EM)[sample(1:(nrow(pbmc)-length(EM)),nrow(pbmc)-length(EM))]),
                  HVGvst = c(hvgvst,setdiff(rownames(pbmc),hvgvst)[sample(1:(nrow(pbmc)-length(hvgvst)),nrow(pbmc)-length(hvgvst))]),
                  HVGdisp = c(hvgdisp,setdiff(rownames(pbmc),hvgdisp)[sample(1:(nrow(pbmc)-length(hvgdisp)),nrow(pbmc)-length(hvgdisp))]),
                  DUBStepR = c(dub,setdiff(rownames(pbmc),dub)[sample(1:(nrow(pbmc)-length(dub)),nrow(pbmc)-length(dub))]),
                  devianceFS = c(devianceFS,setdiff(rownames(pbmc),devianceFS)[sample(1:(nrow(pbmc)-length(devianceFS)),nrow(pbmc)-length(devianceFS))]),
                  trendVar = c(trendvar,setdiff(rownames(pbmc),trendvar)[sample(1:(nrow(pbmc)-length(trendvar)),nrow(pbmc)-length(trendvar))]))
names(pROC_list) <- names(gene.list)
for (i in 1:6){
  pROC_list[[i]] <- roc(c(rep(0,length(g1)),rep(1,length(g2))),
                        # as.numeric(c(g2,g4) %in% gene.list[[i]]),
                        apply(cbind(c(g1,g2)),1,function(x){
                          if (x %in% gene.list[[i]]){
                            which(x==gene.list[[i]])
                          }else{
                            length(gene.list[[i]])+1
                          }
                        }),
                        direction = ">",
                        smoothed = TRUE,
                        # arguments for plot
                        plot=FALSE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=TRUE,
                        print.auc=TRUE, show.thres=TRUE)
}
pdf("./figures/FigureS6_topleft.pdf",width = 6,height = 4)
ggrocs(pROC_list)
dev.off()

# Figure S7 --------------------------------------------------------------
## Figure S7A --------------------------------------------------------------
pbmc <- readRDS("./results/pbmc3k.rds")
ref <- pbmc@active.ident
ref[c("CGGGCATGACCCAA-1","CTTGATTGATCTTC-1")] <- "Platelet"
pbmc <- pbmc[,ref!="Platelet"]
load("./results/pbmc3k_clustering_UMAP.RData")
marker_list <- c("IL7R","CCR7","CD27","SELL","TCF7","CD3E",
                 "CD14", "LYZ",
                 "CD79A","CD37","MS4A1",
                 "GZMK","CD8A","CD8B",
                 "FCGR3A","MS4A7",
                 "NKG7", "GNLY", "CD247","GZMB",
                 "FCER1A", "CST3","CD1C"
                 )
label.tmp <- label.list[[1]]
levels(label.tmp) <- rank(-summary(label.tmp),ties.method = "first")
levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),1:10,c("Memory CD4 T","Naive CD4 T","CD14+ Monocyte",
                                                              "B","CD8 T","FCGR3A+ Monocyte","NK","M-MDSC-like",
                                                              "CD27-CD4+\n Memory T","DC"))
label.tmp <- factor(label.tmp,levels = rev(levels(label.tmp)))
pbmc@active.ident <- label.tmp
pdf("./figures/FigureS7A.pdf",width = 10,height = 6)
DotPlot(pbmc, features = marker_list, cols = c("blue", "red"), 
        dot.scale = 8,idents = levels(label.tmp)[c(-2,-3)]) + RotatedAxis()
dev.off()

## Figure S7B --------------------------------------------------------------
load("./results/pbmc3k_resolution_ARI.RData")
label_ARI_frame <- data.frame()
for (i in 1:nrow(label_ARI)){
  label_ARI_frame <- rbind.data.frame(label_ARI_frame,
                                      data.frame(ARI = label_ARI[i,],
                                                 reso = 0.1*(1:15),
                                                 method = rep(rownames(label_ARI)[i],15)))
}
label_ARI_frame$reso <- factor(label_ARI_frame$reso)
label_ARI_frame$method <- factor(label_ARI_frame$method)
label_ARI_frame$ARI <- round(label_ARI_frame$ARI,2)
pdf("./figures/FigureS7B.pdf",width = 10,height = 4)
ggplot(data = label_ARI_frame,aes(x = reso,y = method,fill = ARI))+
  geom_tile(color = "white",
            lwd = 1,
            linetype = 1)+theme_pubr()+
  scale_y_discrete(limits = c(rownames(label_ARI[-1,])[order(label_ARI[-1,10],decreasing = F)],"Festem"))+
  geom_text(aes(x = reso,y = method,label = ARI),color = "black",size = 5,fontface = "bold")+
  scale_fill_gradient2(low = "white",mid = "#FCD2D1",high = "#FF5C58",midpoint = 0.7) + 
  theme(axis.text.x = element_text(angle = -30),plot.title = element_text(hjust = 0.5),
        panel.border = element_blank()) +
  guides(fill = guide_legend(reverse = F))+
  labs(fill = "Percent",y = NULL,x = "resolution")+
  theme(legend.position = "right")+
  guides(fill = guide_colourbar(label = T,
                                ticks = T,
                                title = NULL))
dev.off()

# Figure S8 ---------------------------------------------------------------
load("./results/pbmc3k_clustering_UMAP.RData")
for (i in 1:length(label.list)){
  label.tmp <- label.list[[i]]
  levels(label.tmp) <- rank(-summary(label.tmp),ties.method = "first")
  if (i == 2) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),3:9,c(4,5,8,6,3,7,10))
  }
  if (i == 3) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),2:8,c(3,2,4,5,6,7,10))
  }
  if (i == 4) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),5:9,c(9,7,6,5,10))
  }
  if (i == 5) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),3:9,c(4,8,5,3,9,7,6))
  }
  if (i == 6) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),2:12,c(3,4,2,6,9,7,5,11,12,10,13))
  }
  umap.tmp <- umap.for.plot(umap.list[[1]],label.tmp)
  
  class_avg <- umap.tmp %>%
    group_by(cluster) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  umap.tmp <- filter(umap.tmp, UMAP_1<0 & UMAP_2>-5)
  umap.tmp$cluster <- factor(umap.tmp$cluster)
  class_avg <- class_avg[class_avg$cluster %in% levels(umap.tmp$cluster),]
  plots.list[[i]] <- ggplot(umap.tmp, aes(x=UMAP_1, y=UMAP_2, color=cluster)) + 
    geom_point(cex=0.1) + theme_pubr()+theme(legend.position="none") +
    # ggrepel::geom_label_repel(data = class_avg,aes(x=UMAP_1,y = UMAP_2,label = allgenes),label.padding = unit(0.1,'cm'))
    geom_text(aes(x=UMAP_1,y = UMAP_2,label = cluster), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 4)+
    labs(title = names(label.list)[i])+
    scale_x_continuous(limits = c(-10,-3))+
    scale_y_continuous(limits = c(-5,8))
  if (i %in% c(2,3)){
    plots.list[[i]] <- plots.list[[i]]+scale_color_manual(values = my.color[c(5)],name = "Clusters")
  } else{
    plots.list[[i]] <- plots.list[[i]]+scale_color_manual(values = my.color[c(9)],name = "Clusters")
  }
}

pbmc <- readRDS("./results/pbmc3k.rds")
ref <- pbmc@active.ident
ref[c("CGGGCATGACCCAA-1","CTTGATTGATCTTC-1")] <- "Platelet"
pbmc <- pbmc[,ref!="Platelet"]
marker_list <- c("CD27","SELL","CCR7","CCL5","CST7","NKG7","CD8A")
feature_plot_list <- vector("list",length = length(marker_list)+2)

for (i in 1:length(marker_list)){
  data.tmp <- plots.list[[1]]$data
  data.tmp$cluster <- factor(data.tmp$cluster)
  exp.tmp <- pbmc@assays$RNA@data[marker_list[i],rownames(data.tmp)]
  data.tmp <- cbind.data.frame(data.tmp,exp = exp.tmp)
  umap.tmp$cluster <- factor(umap.tmp$cluster)
  feature_plot_list[[i]] <- 
    ggplot(data.tmp,mapping = aes(x=UMAP_1, y=UMAP_2,color = exp)) + 
    geom_point(cex=0.5,shape = 16)+ 
    scale_color_gradient(low = "grey75",high = "#CD0404")+
    theme_pubr()+theme(legend.position="right") +
    scale_x_continuous(limits = c(-10,-3))+
    scale_y_continuous(limits = c(-5,8))+
    labs(title = marker_list[i])
}
feature_plot_list[8:9] <- plots.list[c(1,5)]

pdf("./figures/FigureS8.pdf",height = 8,width = 8)
ggarrange(plotlist = feature_plot_list,ncol = 3,nrow = 3)
dev.off()

# Figure S18 --------------------------------------------------------------
load("./results/pbmc3k_clustering_UMAP.RData")
pbmc <- readRDS("./results/pbmc3k.rds")
ref <- pbmc@active.ident
ref[c("CGGGCATGACCCAA-1","CTTGATTGATCTTC-1")] <- "Platelet"
pbmc <- pbmc[,ref!="Platelet"]
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes,match = rownames(pbmc))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes,match = rownames(pbmc))
pbmc <- CellCycleScoring(pbmc,g2m.features = g2m_genes,
                         s.features = s_genes)

data.tmp <- plots.list[[1]]$data
data.tmp$cluster <- factor(data.tmp$cluster)
exp.tmp <- pbmc@assays$RNA@data["MAP1LC3B",]
data.tmp <- cbind.data.frame(data.tmp,exp = exp.tmp)
p1 <- ggplot(mapping = aes(x=UMAP_1, y=UMAP_2)) + 
  geom_point(data = data.tmp,
             aes(color=exp),
             cex=0.2,shape = 16)+ 
  scale_color_gradient(low = "grey75",high = "blue")+
  theme_pubr()+theme(legend.position="none") +
  labs(title = "MAP1LC3B")
marker_exp <- data.frame(exp = as.numeric(pbmc@assays$RNA@data["MAP1LC3B",]),
                         id = pbmc@meta.data$Phase)
set.seed(321)
noise <- rnorm(n = nrow(marker_exp)) / 100000
marker_exp$exp <- marker_exp$exp+noise
color <- hue_pal()(3)
p2 <- ggplot(marker_exp, aes(x = factor(id,levels = c("G1","G2M","S")),
                             y = exp, fill = id)) +
  geom_violin(scale = "width", adjust = 0.8, trim = T) +
  scale_y_continuous(expand = c(0, 0), position="left", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
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
  ggtitle("MAP1LC3B")+xlab("") + ylab("Expression Level")
pdf("FigureS18_MAP1LC3B.pdf",width = 6,height = 3)
ggarrange(p1,p2,nrow = 1,labels = c("A","B"))
dev.off()

data.tmp <- plots.list[[1]]$data
data.tmp$cluster <- factor(data.tmp$cluster)
exp.tmp <- pbmc@assays$RNA@data["HMGB2",]
data.tmp <- cbind.data.frame(data.tmp,exp = exp.tmp)
p1 <- ggplot(mapping = aes(x=UMAP_1, y=UMAP_2)) + 
  geom_point(data = data.tmp,
             aes(color=exp),
             cex=0.2,shape = 16)+ 
  scale_color_gradient(low = "grey75",high = "blue")+
  theme_pubr()+theme(legend.position="none") +
  labs(title = "HMGB2")
marker_exp <- data.frame(exp = as.numeric(pbmc@assays$RNA@data["HMGB2",]),
                         id = pbmc@meta.data$Phase)
set.seed(321)
noise <- rnorm(n = nrow(marker_exp)) / 100000
marker_exp$exp <- marker_exp$exp+noise
color <- hue_pal()(3)
p2 <- ggplot(marker_exp, aes(x = factor(id,levels = c("G1","G2M","S")),
                             y = exp, fill = id)) +
  geom_violin(scale = "width", adjust = 0.8, trim = T) +
  scale_y_continuous(expand = c(0, 0), position="left", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
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
  ggtitle("HMGB2")+xlab("") + ylab("Expression Level")
pdf("FigureS18_HMGB2.pdf",width = 6,height = 3)
ggarrange(p1,p2,nrow = 1,labels = c("A","B"))
dev.off()