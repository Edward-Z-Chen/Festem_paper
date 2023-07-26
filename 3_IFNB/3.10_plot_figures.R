library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(getopt)
my.color <- hue_pal()(17)
names(my.color) <- 1:17
tsne.for.plot <- function(tsne,cluster){
  tsne <- tsne@cell.embeddings
  tsne <- as.data.frame(tsne)
  cbind(tsne,cluster = cluster)
}
if (!file.exists("./figures")){
  dir.create("./figures")
}

spec <- matrix(
  c("all_analysis", "a", 0,"logical", "Perform all analysis"
  ),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)

# Figure 3 (A) -- top right and Figure S5 second from top ------------------------------------------------
if (!is.null(opt$all_analysis)){
## Summarizing results from various DEG detection methods
ifnb <- readRDS("./results/ifnb_ctrl.rds")
results <- matrix(nrow = nrow(ifnb),ncol = 19,
                  dimnames = list(rownames(ifnb),
                                  c("Festem","Festem_gamma0.01","Festem_13g","Festem_17g",
                                    "DEseq2", "DEseq2-full",       
                                    "EdgeR", "EdgeR-full","MAST",              
                                    "Wilcoxon","MAST-f",            
                                    "Wilcoxon-f","FC",
                                    "singleCellHaystack-PCA",
                                    "singleCellHaystack-TSNE",
                                    "singleCellHaystack-PCA-20pc",
                                    "singleCellHaystack-TSNE-20pc","TN_test","ROSeq")))
FC.list <- vector("list",nlevels(ifnb@meta.data$seurat_annotations))
for (i in 1:length(FC.list)){
  FC.list[[i]] <- FoldChange(ifnb,ident.1 = levels(ifnb@meta.data$seurat_annotations)[i],
                             group.by = "seurat_annotations")[,1]
}
FC.list <- matrix(unlist(FC.list),ncol = nlevels(ifnb@meta.data$seurat_annotations),byrow = F)
FC <- apply(FC.list,1,function(x){max(abs(x))})
results[,"FC"] <- FC
load("./results/ifnb_ctrl_Festem.RData")
results[colnames(em.result),"Festem"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result.9)[em.result.9[2,]<0],"Festem"] <- NA
load("./results/ifnb_ctrl_Festem_gamma0.01.RData")
results[colnames(em.result),"Festem_gamma0.01"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result.9)[em.result.9[2,]<0],"Festem_gamma0.01"] <- NA
load("./results/ifnb_ctrl_Festem_13g.RData")
results[colnames(em.result),"Festem_13g"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result.9)[em.result.9[2,]<0],"Festem_13g"] <- NA
load("./results/ifnb_ctrl_Festem_17g.RData")
results[colnames(em.result),"Festem_17g"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result.9)[em.result.9[2,]<0],"Festem_17g"] <- NA
load("./results/ifnb_DEseq2.RData")
results[DEseq.results@rownames,"DEseq2"] <- DEseq.results@listData$padj
results[DEseq.results@rownames,"DEseq2-full"] <- p.adjust(DEseq.results@listData$pvalue,"BH")
load("./results/ifnb_edgeR.RData")
results[names(EdgeR.result),"EdgeR"] <- p.adjust(EdgeR.result,"BH")
load("./results/ifnb_edgeR_full.RData")
results[names(EdgeR.result),"EdgeR-full"] <- p.adjust(EdgeR.result,"BH")
load("./results/ifnb_mast.RData")
results[rownames(mast.results),"MAST"] <- p.adjust(mast.results[,3],"BH")
load("./results/ifnb_wilcoxon.RData")
wil.result <- matrix(unlist(wil.result),ncol = length(wil.result),byrow = F)
wil.result <- apply(wil.result,1,function(x){min(x)*length(wil.result)})
results[,"Wilcoxon"] <- p.adjust(wil.result,"BH")
results[,"Wilcoxon-f"] <- results[,"Wilcoxon"]
results[results[,"FC"]<=0.2,"Wilcoxon-f"] <- NA
results[,"MAST-f"] <- results[,"MAST"]
results[results[,"FC"]<=0.2,"MAST-f"] <- NA
load("./results/ifnb_haystack.RData")
results[,"singleCellHaystack-PCA"] <- exp(haystack_pca$results$log.p.adj)
results[,"singleCellHaystack-TSNE"] <- exp(haystack_umap$results$log.p.adj)
load("./results/ifnb_haystack_20PC.RData")
results[,"singleCellHaystack-PCA-20pc"] <- exp(haystack_pca$results$log.p.adj)
results[,"singleCellHaystack-TSNE-20pc"] <- exp(haystack_umap$results$log.p.adj)
load("./results/ifnb_ROSeq.RData")
results[,"ROSeq"] <- p.adjust(roseq.tmp,"BH")

tn_test <- read.csv("./results/ifnb_TN_test.csv",
                    header = F)
tn_test <- as.matrix(tn_test)
tn_test <- apply(tn_test,1,function(x){min(x)*length(x)})
results[,"TN_test"] <- p.adjust(tn_test,"BH")
save(results,file = "./results/ifnb_ctrl_DEG_results.RData")
} else{
  load("./results/ifnb_ctrl_DEG_results.RData")
}
## Figure S5 second from top (left) -----------------------------------------------------------
calc_house_percent <- function(DEG_list,house_list){
  sum(house_list%in%DEG_list)/length(DEG_list)
}
calc_house_retrieve <- function(DEG_list,house_list){
  sum(DEG_list%in%house_list)/length(house_list)
}
load("./results/ifnb_ctrl_silver_standard.RData")
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
names(my_color) <- c("Festem","Festem_gamma0.01","Festem_13g","Festem_17g",
                     "DEseq2","DEseq2-full","EdgeR","EdgeR-full",
                     "MAST","Wilcoxon","MAST-f","Wilcoxon-f","ROSeq",
                     "singleCellHaystack-PCA","singleCellHaystack-TSNE",
                     "singleCellHaystack-PCA-20pc","singleCellHaystack-TSNE-20pc","TN_test")
my_shape <- c(16,15,17,18,16,7,16,7,16,16,9,9,16,10,16,12,13,16)
names(my_shape) <- names(my_color)

pdf("./figures/FigureS5_ifnb_left.pdf",height = 6,width = 6)
ggplot(data = power_frame2,mapping = aes(x = recall, y = precision))+
  geom_point(aes(color = method, shape = method),
             size = 4,stroke = 1)+
  scale_x_continuous(limits = c(0.2,1))+
  scale_y_continuous(limits = c(0.1,1))+
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

## Figure S5 second from top (right) -----------------------------------------------------------
num_reject <- apply(p_values_full,2,function(x){sum(x<0.05,na.rm = T)})
num_reject <- data.frame(method = colnames(p_values_full),
                         num = num_reject)
num_reject <- num_reject[order(num_reject$num,decreasing = F),]
num_reject$method <- factor(num_reject$method,levels = num_reject$method)

pdf("./figures/FigureS5_ifnb_right.pdf",height = 8,width = 8)
ggplot(num_reject,aes(x = method,y = num,fill = method))+
  geom_bar(stat = "identity")+
  theme_pubr()+
  geom_text(aes(label=num), vjust=1.6, color="white", size=3.5)+
  scale_fill_manual(values = my_color)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90))+
  guides(color = guide_legend(nrow = 4))
dev.off()

## Figure 3 (A) -- top right ------------------------------------------------
power_frame2 <- power_frame2[c(1,5,7,9,10,14,17,18),]
power_frame2[c(1,6),1] <- c("Festem","singleCellHaystack")
my_color <- my_color[c(1,5,7,9,10,13,15,18)]
names(my_color)[7] <- "singleCellHaystack"
my_shape <- c(16,15,17,18,7,8,9,10)
names(my_shape) <- names(my_color)

pdf("./figures/Figure3A_topright.pdf",height = 4,width = 3.5)
ggplot(data = power_frame2,mapping = aes(x = recall, y = precision))+
  geom_point(aes(color = method, shape = method),
             size = 4,stroke = 1)+
  scale_x_continuous(limits = c(0.1,1))+
  scale_y_continuous(limits = c(0.1,1))+
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

num_reject <- num_reject[num_reject$method %in% c(power_frame2$method,"singleCellHaystack-TSNE"),]
num_reject$method <- factor(num_reject$method)
levels(num_reject$method) <- plyr::mapvalues(levels(num_reject$method),
                                             c("singleCellHaystack-TSNE"),
                                             c("singleCellHaystack"))

pdf("./figures/Figure3A_topright_inset.pdf",height = 4,width = 4)
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

# Figure 4 & Figure S11B ----------------------------------------------------------------
load("./results/ifnb_ctrl_clustering_tSNE.RData")
my.color <- hue_pal()(17)
names(my.color) <- 1:17

## Figure 4E ----------------------------------------------------------------
label.tmp <- label.list[[1]]
levels(label.tmp) <- rank(-summary(label.tmp),ties.method = "first")
tsne.tmp <- tsne.for.plot(tsne.list[[1]],label.tmp)

class_avg <- tsne.tmp %>%
  group_by(cluster) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )
class_avg <- cbind(class_avg,class_avg[,1])
colnames(class_avg)[4] <- "cluster.anno"
class_avg[,4] <- plyr::mapvalues(class_avg[,4],1:17,c("CD14 Monocyte","CD4 Memory T","CD4 Naive T",
                                                      "CD14 Monocyte","CD16 monocyte","CD14 Monocyte",
                                                      "B","T activated","CD8 T","NK","DC", "B Activated",
                                                      "T cell:monocyte complex","HSP+ CD4 T",
                                                      "IFNhi CD14 Monocyte","	CD34+ progenitors","pDC"))
levels(tsne.tmp$cluster)[c(1,4,6)] <- 1
pdf("./figures/Figure4E.pdf",width = 4,height = 4)
ggplot(tsne.tmp, aes(x=tSNE_1, y=tSNE_2, color=cluster)) + 
  geom_point(cex=0.1) + theme_pubr()+
  theme(legend.position="none",text = element_text(size = 15)) +
  # ggrepel::geom_label_repel(data = class_avg,aes(x=tSNE_1,y = tSNE_2,label = cluster.anno),label.padding = unit(0.1,'cm'),label.size = 0.15)+
  geom_text(aes(x=tSNE_1,y = tSNE_2,label = cluster.anno), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 3)+
  labs(title = names(label.list)[1])+scale_color_manual(values = my.color,name = "Clusters")
dev.off()

## Figure 4F ----------------------------------------------------------------
library(cowplot)
ifnb <- readRDS("./results/ifnb_ctrl.rds")
ifnb <- NormalizeData(ifnb)
ifnb <- ifnb[,!ifnb@meta.data$seurat_annotations%in%c("Eryth","Mk")]
marker_list <- c("ISG15","ISG20","IFI6","IFIT1","IFITM3","IFIT3","IFI35",
                 "LY6E","APOBEC3A","RSAD2","OAS1","MX1","TNFSF10","TNFSF13B")
ifnb <- ScaleData(ifnb,features = marker_list)
marker_exp <- data.frame(exp = as.numeric(ifnb@assays$RNA@scale.data[marker_list,plots.list[[1]]$data$cluster=="14"]),
                         gene = rep(marker_list,sum(plots.list[[1]]$data$cluster=="14")),
                         id = "IFNhi CD14 Monocyte")
marker_exp <- rbind(marker_exp,data.frame(exp = as.numeric(ifnb@assays$RNA@scale.data[marker_list,plots.list[[1]]$data$cluster %in% c(0,3,5)]),
                                          gene = rep(marker_list,sum(plots.list[[1]]$data$cluster %in% c(0,3,5))),
                                          id = "CD14 Monocyte"))
set.seed(321)
noise <- rnorm(n = nrow(marker_exp)) / 100000
marker_exp$exp <- marker_exp$exp+noise
marker_exp$gene <- factor(marker_exp$gene,levels = marker_list)
color <- my.color[c(1,15)]
names(color) <- c("CD14 Monocyte","IFNhi CD14 Monocyte")

pdf("./figures/Figure4F.pdf",width = 4,height = 5.5)
ggplot(marker_exp, aes(x = factor(id,levels = c("IFNhi CD14 Monocyte","CD14 Monocyte")),
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
  ggtitle("Festem")+xlab("") + ylab("Expression Level")
dev.off()

## Figure 4G --------------------------------------------------------------
library(cowplot)
ifnb <- readRDS("./results/ifnb_ctrl.rds")
ifnb <- NormalizeData(ifnb)
ifnb <- ifnb[,!ifnb@meta.data$seurat_annotations%in%c("Eryth","Mk")]
marker_list <- c("HSPA1A" ,"HSPB1" ,"HSPA1B" ,"HSPH1" ,
                 "HSPE1" ,"HSPD1" ,"HSP90AB1" ,"DNAJB1",
                 "HSP90AA1" ,"DNAJA1","HSPA8")
ifnb <- ScaleData(ifnb,features = marker_list)
marker_exp <- data.frame(exp = as.numeric(ifnb@assays$RNA@scale.data[marker_list,plots.list[[1]]$data$cluster=="3"]),
                         gene = rep(marker_list,sum(plots.list[[1]]$data$cluster=="3")),
                         id = "CD4 Naive T")
marker_exp <- rbind(marker_exp,data.frame(exp = as.numeric(ifnb@assays$RNA@scale.data[marker_list,plots.list[[1]]$data$cluster=="2"]),
                                          gene = rep(marker_list,sum(plots.list[[1]]$data$cluster=="2")),
                                          id = "CD4 Memory T"))
marker_exp <- rbind(marker_exp,data.frame(exp = as.numeric(ifnb@assays$RNA@scale.data[marker_list,plots.list[[1]]$data$cluster=="9"]),
                                          gene = rep(marker_list,sum(plots.list[[1]]$data$cluster=="9")),
                                          id = "CD8 T"))
marker_exp <- rbind(marker_exp,data.frame(exp = as.numeric(ifnb@assays$RNA@scale.data[marker_list,plots.list[[1]]$data$cluster=="8"]),
                                          gene = rep(marker_list,sum(plots.list[[1]]$data$cluster=="8")),
                                          id = "T Activated"))
marker_exp <- rbind(marker_exp,data.frame(exp = as.numeric(ifnb@assays$RNA@scale.data[marker_list,plots.list[[1]]$data$cluster=="14"]),
                                          gene = rep(marker_list,sum(plots.list[[1]]$data$cluster=="14")),
                                          id = "HSP+ CD4 T"))
set.seed(321)
noise <- rnorm(n = nrow(marker_exp)) / 100000
marker_exp$exp <- marker_exp$exp+noise
marker_exp$gene <- factor(marker_exp$gene,levels = marker_list)

my.color <- hue_pal()(17)
names(my.color) <- 1:17
color <- my.color[c(3,2,9,8,14)]
names(color) <- unique(marker_exp$id)

pdf("./figures/Figure4G.pdf",width = 4,height = 4)
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
  ggtitle("HSP+ CD4 T (Cluster 12)")+xlab("") + ylab("Expression Level")
dev.off()

## Figure S11B ----------------------------------------------------------------
ifnb <- readRDS("./results/ifnb_ctrl.rds")
ifnb <- NormalizeData(ifnb)
ifnb <- ifnb[,!ifnb@meta.data$seurat_annotations%in%c("Eryth","Mk")]
marker_list <- list(c("ISG15","ISG20","IFI6","IFIT1","IFITM3","IFIT3","IFI35",
                      "LY6E","APOBEC3A","RSAD2","OAS1","MX1","TNFSF10","TNFSF13B"))
ifnb <- AddModuleScore(ifnb,marker_list,name = "IFN")

load("./results/ifnb_ctrl_clustering_tSNE.RData")
data.tmp <- plots.list[[1]]$data
exp.tmp <- ifnb@meta.data[,"IFN1"]
data.tmp <- cbind.data.frame(data.tmp,exp = exp.tmp)

p_ctrl <- ggplot(mapping = aes(x=tSNE_1, y=tSNE_2)) + 
  geom_point(data = data.tmp,
             aes(color=exp),
             cex=0.1,shape = 16)+ 
  scale_color_gradient(low = "grey75",high = "#CD0404")+
  theme_pubr()+theme(legend.position="right") +
  labs(title = "Control",color = "IFN score")

load("./results/ifnb_stim_clustering_tSNE.RData")
ifnb_stim <- readRDS("./results/ifnb_stim.rds")
ifnb_stim <- NormalizeData(ifnb_stim)
ifnb_stim <- ifnb_stim[,!ifnb_stim@meta.data$seurat_annotations%in%c("Eryth","Mk")]
ifnb_stim <- AddModuleScore(ifnb_stim,marker_list,name = "IFN")
data.tmp <- tsne.for.plot(tsne,label)
exp.tmp <- ifnb_stim@meta.data[,"IFN1"]
data.tmp <- cbind.data.frame(data.tmp,exp = exp.tmp)
p_stim <- 
  ggplot(mapping = aes(x=tSNE_1, y=tSNE_2)) + 
  geom_point(data = data.tmp,
             aes(color=exp),
             cex=0.1,shape = 16)+ 
  scale_color_gradient(low = "grey75",high = "#CD0404")+
  theme_pubr()+theme(legend.position="right") +
  labs(title = "Stimulated",color = "IFN score")


pdf("./figures/S11B.pdf",width = 7,height = 4)
ggarrange(p_ctrl,p_stim,ncol = 2, labels = c("A","B"),common.legend = T,legend = "bottom")
dev.off()

# Figure S7 -- top right ---------------------------------------------------
load("./results/ifnb_ctrl_hvggenes.RData")
ifnb <- readRDS("./results/ifnb_ctrl.rds")
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb,nfeatures = nrow(ifnb))
hvgvst <- VariableFeatures(ifnb)
ifnb <- FindVariableFeatures(ifnb,nfeatures = nrow(ifnb),selection.method = "disp")
hvgdisp <- VariableFeatures(ifnb)

load("./results/ifnb_ctrl_silver_standard.RData")
g1 <- rownames(moran_h)[abs(moran_h[,1])<0.02 & moran_h[,2]>0.05]
moran_nonh <- moran_nonh[!is.na(moran_nonh[,1]),]
g2 <- rownames(moran_nonh)[moran_nonh[,1]>0.1]
g2 <- g2[gene_percent[g2]>0.05]
load("./results/ifnb_ctrl_clustering_tSNE.RData")
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
gene.list <- list(Festem = c(EM,setdiff(rownames(ifnb),EM)[sample(1:(nrow(ifnb)-length(EM)),nrow(ifnb)-length(EM))]),
                  HVGvst = c(hvgvst,setdiff(rownames(ifnb),hvgvst)[sample(1:(nrow(ifnb)-length(hvgvst)),nrow(ifnb)-length(hvgvst))]),
                  HVGdisp = c(hvgdisp,setdiff(rownames(ifnb),hvgdisp)[sample(1:(nrow(ifnb)-length(hvgdisp)),nrow(ifnb)-length(hvgdisp))]),
                  DUBStepR = c(dub,setdiff(rownames(ifnb),dub)[sample(1:(nrow(ifnb)-length(dub)),nrow(ifnb)-length(dub))]),
                  devianceFS = c(devianceFS,setdiff(rownames(ifnb),devianceFS)[sample(1:(nrow(ifnb)-length(devianceFS)),nrow(ifnb)-length(devianceFS))]),
                  trendVar = c(trendvar,setdiff(rownames(ifnb),trendvar)[sample(1:(nrow(ifnb)-length(trendvar)),nrow(ifnb)-length(trendvar))]))
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
pdf("./figures/FigureS7_topright.pdf",width = 6,height = 4)
ggrocs(pROC_list)
dev.off()

# Figure S10 ---------------------------------------------------------------
load("./results/ifnb_ctrl_clustering_tSNE.RData")
ifnb <- readRDS("./results/ifnb_ctrl.rds")
ifnb <- NormalizeData(ifnb)
ifnb <- ifnb[,!ifnb@meta.data$seurat_annotations%in%c("Eryth","Mk")]
marker_list <- c("CD3D", "CREM", "IL7R","CCR7","CD27",
                 "SELL", "GIMAP5", "CACYBP", "TCF7",
                 "GNLY", "NKG7", "CCL5", "CD247","GZMB",
                 "CD8A", 
                 "MS4A1", "CD79A", "CD37",
                 "MIR155HG", "NME1", 
                 "FCGR3A", "VMO1", "MS4A7",
                 "CCL2", "S100A9", "CD14", "LYZ",
                 "HLA-DQA1", "GPR183", "FCER1A","CST3","CD1C",
                 "TSPAN13", "IL3RA", "IGJ",
                 "HSPA1A","HSPB1","HSPA1B","HSPH1","HSPE1","HSPD1","HSPA6",
                 "CCL7","CCL3","CCL4",
                 "CD34","TPSAB1","GATA2","SNHG7")
label.tmp <- label.list[[1]]
levels(label.tmp) <- rank(-summary(label.tmp),ties.method = "first")
levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),1:17,c("CD14 Monocyte","CD4 Memory T","CD4 Naive T",
                                                              "CD14 Monocyte","CD16 monocyte","CD14 Monocyte",
                                                              "B","T activated","CD8 T","NK","DC", "B Activated",
                                                              "T cell:monocyte complex","HSP+ CD4 T",
                                                              "IFNhi CD14 Monocyte","CD34+ progenitors","pDC"))
label.tmp <- factor(label.tmp,levels = c("T cell:monocyte complex",
                                         "IFNhi CD14 Monocyte","CD34+ progenitors",
                                         "HSP+ CD4 T",
                                         "pDC",
                                         "DC","CD14 Monocyte","CD16 monocyte","B Activated",
                                         "B","CD8 T","NK","T activated",
                                         "CD4 Naive T",
                                         "CD4 Memory T"
                                         ))
ifnb@active.ident <- label.tmp
ifnb <- ScaleData(ifnb,features = marker_list)
pdf("./figures/FigureS10.pdf",width = 15,height = 5)
DotPlot(ifnb, features = marker_list, cols = c("blue", "red"), 
        dot.scale = 8,idents = levels(label.tmp)) + RotatedAxis()
dev.off()

# Figure S11 --------------------------------------------------------------
load("./results/ifnb_ctrl_clustering_tSNE.RData")
set.seed(321)
my.color <- hue_pal()(35)
names(my.color) <- sample(1:35,size = 35)
for (i in 1:length(label.list)){
  label.tmp <- label.list[[i]]
  levels(label.tmp) <- rank(-summary(label.tmp),ties.method = "first")
  if (i == 2) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),1:18,
                                         c(4,2,1,18,5,3,6,7,9,8,10,11,12,15,19,20,17,16))
  }
  if (i == 3) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),1:16,c(1,2,4,3,5,8,7,9,10,11,13,12,15,20,17,16))
  }
  if (i == 5) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),1:18,c(2,1,3,5,4,6,7,8,9,21,10,11,12,13,14,16,22,17))
  }
  if (i == 6) {
    levels(label.tmp) <- plyr::mapvalues(levels(label.tmp),1:16,c(3,2,1,4,5,7,8,6,9,11,10,13,12,20,15,17))
  }
  tsne.tmp <- tsne.for.plot(tsne.list[[1]],label.tmp)
  
  class_avg <- tsne.tmp %>%
    group_by(cluster) %>%
    summarise(
      tSNE_1 = median(tSNE_1),
      tSNE_2 = median(tSNE_2)
    )
  plots.list[[i]] <- ggplot(tsne.tmp, aes(x=tSNE_1, y=tSNE_2, color=cluster)) + 
    geom_point(cex=0.1) + theme_pubr()+theme(legend.position="none",text = element_text(size = 15)) +
    # ggrepel::geom_label_repel(data = class_avg,aes(x=tSNE_1,y = tSNE_2,label = allgenes),label.padding = unit(0.1,'cm'))
    geom_text(aes(x=tSNE_1,y = tSNE_2,label = cluster), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 4)+
    labs(title = names(label.list)[i])+scale_color_manual(values = my.color,name = "Clusters")
}

pdf("./figures/FigureS11.pdf",width = 10, height = 7,onefile = T)
ggarrange(plotlist = plots.list,nrow = 2,ncol = 3,legend = "none")
dev.off()

# ## Figure S10B --------------------------------------------------------------
# load("./results/ifnb_ctrl_resolution_ARI.RData")
# label_ARI_frame <- data.frame()
# for (i in 1:nrow(label_ARI)){
#   label_ARI_frame <- rbind.data.frame(label_ARI_frame,
#                                       data.frame(ARI = label_ARI[i,],
#                                                  reso = 0.1*(1:15)+0.5,
#                                                  method = rep(rownames(label_ARI)[i],15)))
# }
# label_ARI_frame$reso <- factor(label_ARI_frame$reso)
# label_ARI_frame$method <- factor(label_ARI_frame$method)
# label_ARI_frame$ARI <- round(label_ARI_frame$ARI,2)
# 
# pdf("./figures/FigureS10B.pdf",width = 8,height = 4)
# ggplot(data = label_ARI_frame,aes(x = reso,y = method,fill = ARI))+
#   geom_tile(color = "white",
#             lwd = 1,
#             linetype = 1)+theme_pubr()+
#   scale_y_discrete(limits = c(rownames(label_ARI[-1,])[order(label_ARI[-1,10],decreasing = F)],"Festem"))+
#   geom_text(aes(x = reso,y = method,label = ARI),color = "black",size = 5,fontface = "bold")+
#   scale_fill_gradient2(low = "white",mid = "#FCD2D1",high = "#FF5C58",midpoint = 0.5) + 
#   theme(axis.text.x = element_text(angle = -30),plot.title = element_text(hjust = 0.5),
#         panel.border = element_blank()) +
#   guides(fill = guide_legend(reverse = F))+
#   labs(fill = "Percent",y = NULL,x = "resolution")+
#   theme(legend.position = "right")+
#   guides(fill = guide_colourbar(label = T,
#                                 ticks = T,
#                                 title = NULL))
# dev.off()


# Figure S12A --------------------------------------------------------------
load("./results/ifnb_stim_clustering_tSNE.RData")
ifnb_stim <- readRDS("./results/ifnb_stim.rds")
ifnb_stim <- NormalizeData(ifnb_stim)
ifnb_stim <- ifnb_stim[,!ifnb_stim@meta.data$seurat_annotations%in%c("Eryth","Mk")]

label.tmp <- label
levels(label.tmp) <- rank(-summary(label.tmp),ties.method = "first")
tsne.tmp <- tsne.for.plot(tsne,label.tmp)

class_avg <- tsne.tmp %>%
  group_by(cluster) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )
class_avg <- cbind(class_avg,class_avg[,1])
colnames(class_avg)[4] <- "cluster.anno"
class_avg[,4] <- plyr::mapvalues(class_avg[,4],1:16,c("IFNhi CD14 Monocyte","CD4 Naive T","CD4 Memory T","B",
                                                      "IFNhi CD14 Monocyte","CD16 Monocyte","CD8 T",
                                                      "IFNhi CD14 Monocyte","T activated","CD4 Naive T","NK",
                                                      "DC","B Activated","T cell:monocyte complex","T activated","pDC"
                                                      ))
levels(tsne.tmp$cluster) <- plyr::mapvalues(levels(tsne.tmp$cluster),c(1,5,8,2,10,9,15),c(1,1,1,10,10,15,15))

pdf("./figures/FigureS12A.pdf",width = 4,height = 4)
ggplot(tsne.tmp, aes(x=tSNE_1, y=tSNE_2, color=cluster)) + 
  geom_point(cex=0.1) + theme_pubr()+
  theme(legend.position="none",text = element_text(size = 15)) +
  # ggrepel::geom_label_repel(data = class_avg,aes(x=tSNE_1,y = tSNE_2,label = cluster.anno),label.padding = unit(0.1,'cm'),label.size = 0.15)+
  geom_text(aes(x=tSNE_1,y = tSNE_2,label = cluster.anno), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 3)+
  labs(title = "Stimulated")+scale_color_manual(values = my.color,name = "Clusters")
dev.off()

# Figure S13 --------------------------------------------------------------
ifnb <- SeuratData::LoadData("ifnb")
ifnb <- NormalizeData(ifnb)
ifnb <- ifnb[,!ifnb@meta.data$seurat_annotations%in%c("Eryth","Mk")]
tmp <- data.frame(gene = rownames(ifnb),
                  ctrl_p = NA, ctrl_EM = NA, stim_p = NA,stim_EM = NA)
rownames(tmp) <- rownames(ifnb)
load("./results/ifnb_ctrl_Festem.RData")
tmp[colnames(em.result),2] <- em.result[1,]
tmp[colnames(em.result.9),3] <- em.result.9[2,]
load("./results/ifnb_stim_Festem.RData")
tmp[colnames(em.result),4] <- em.result[1,]
tmp[colnames(em.result.9),5] <- em.result.9[2,]

result.EM = data.frame(names = tmp$gene, 
                       p = NA, 
                       EM = NA)
for (i in 1:nrow(result.EM)){
  if (sum(is.na(tmp[i,c(2,4)])) < 2){
    result.EM[i,2] <- min(tmp[i,2],tmp[i,4]) * 2
    result.EM[i,3] <- max(tmp[i,3],tmp[i,5])
  }
}
result.EM$p <- p.adjust(result.EM$p,"BH")
tmp.na <- result.EM[is.na(result.EM$p),1]
result.EM <- result.EM[!is.na(result.EM$p),]
gene.names <- result.EM[result.EM$p<0.05 & result.EM$EM>0,]
gene.names <- gene.names[order(-gene.names$EM),]
tmp <- result.EM[result.EM$p>=0.05 & result.EM$EM>0,]
tmp <- tmp[order(-tmp$EM),]
gene.names <- rbind(gene.names,tmp)
tmp <- result.EM[result.EM$EM<=0,]
tmp <- tmp[order(tmp$p,-tmp$EM),]
gene.names <- rbind(gene.names,tmp)
gene.names <- gene.names[,1]
gene.names <- c(gene.names,tmp.na)
EM <- gene.names

ifnb <- ScaleData(ifnb,features = EM[1:2500])
ifnb <- RunPCA(ifnb, verbose = FALSE,features = EM[1:2500])
ifnb <- harmony::RunHarmony(ifnb,"stim",plot_convergence = T, lambda = 30)
ifnb <- FindNeighbors(object = ifnb, dims = 1:20,reduction = "harmony")
ifnb <- RunTSNE(ifnb,dims = 1:20,reduction = "harmony")
cell_type <- rep(NA,ncol(ifnb))
names(cell_type) <- colnames(ifnb)
load("./results/ifnb_ctrl_clustering_tSNE.RData")
levels(label.list[[1]]) <- plyr::mapvalues(levels(label.list[[1]]),c(0,5,3),c(0,0,0))
cell_type[names(label.list[[1]])] <- paste0("C_",label.list[[1]])
load("./results/ifnb_stim_clustering_tSNE.RData")
levels(label) <- plyr::mapvalues(levels(label),c(1,5,8,2,10,9,15)-1,c(1,1,1,10,10,15,15)-1)
cell_type[names(label)] <- paste0("S_",label)
ifnb <- AddMetaData(ifnb,cell_type,"anno")
DimPlot(ifnb,label = T,group.by = "stim")
DimPlot(ifnb,label = T,group.by = "anno")

tsne.tmp <- tsne.for.plot(ifnb@reductions[["tsne"]],ifnb@meta.data$anno)
class_avg <- tsne.tmp %>%
  group_by(cluster) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )
p1 <- ggplot(tsne.tmp, aes(x=tSNE_1, y=tSNE_2, color=cluster)) + 
  geom_point(cex=0.1) + theme_pubr()+theme(legend.position="right")
p2 <- ggplot(tsne.tmp, aes(x=tSNE_1, y=tSNE_2, color=cluster)) + 
  geom_point(cex=0.1) + theme_pubr()+theme(legend.position="right")+
  scale_color_manual(values = c("C_14" = "red", "S_0" = "lightblue","C_0"="lightgreen"))

pdf("./figures/FigureS13.pdf",width = 12.5,height = 5)
ggarrange(p1,p2,ncol = 2,widths = c(1.1,1))
dev.off()

# ## Figure S12C --------------------------------------------------------------
# load("./results/ifnb_ctrl_clustering_tSNE.RData")
# 
# pdf("./figures/FigureS12C.pdf",width = 6,height = 4)
# list(Festem = gene.list[[1]],devianceFS = gene.list[[5]]) %>%
#   ggVennDiagram::ggVennDiagram(label_alpha = 0)+
#   scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
# dev.off()
# 
# ## Figure S12D -------------------------------------------------------------
# load("./results/ifnb_ctrl_silver_standard.RData")
# moran <- rbind(moran_h,moran_nonh)
# x1 <- moran[setdiff(gene.list[[1]],gene.list[[5]]),]
# x2 <- moran[setdiff(gene.list[[5]],gene.list[[1]]),]
# 
# pdf("./figures/FigureS12D.pdf",width = 4,height = 4)
# data.frame(name = factor(c(rep("Festem only",nrow(x1)),rep("devianceFS only",nrow(x2))),
#                          levels = c("Festem only","devianceFS only")),
#            value = c(x1$observed,x2$observed)) %>%
#   ggplot( aes(x=name, y=value, fill=name)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("Festem only" = "red","devianceFS only" = "#FF7F00"))+
#   theme_pubr() +
#   theme(
#     legend.position="none",
#     plot.title = element_text(size=11)
#   ) +
#   stat_compare_means()+
#   xlab("")+ylab("Moran index")
# dev.off()
# 
# ## Figure S12E -------------------------------------------------------------
# load("./results/ifnb_ctrl_silver_standard.RData")
# moran <- rbind(moran_h,moran_nonh)
# x1 <- moran[setdiff(gene.list[[1]],gene.list[[2]]),]
# x2 <- moran[setdiff(gene.list[[2]],gene.list[[1]]),]
# 
# p1 <- data.frame(name = c(rep("Festem only",nrow(x1)),rep("HVGvst only",nrow(x2))),
#            value = c(x1$observed,x2$observed)) %>%
#   ggplot( aes(x=name, y=value, fill=name)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("Festem only" = "red","HVGvst only" = "#6A3D9A"))+
#   theme_pubr() +
#   theme(
#     legend.position="none",
#     plot.title = element_text(size=11)
#   ) +
#   stat_compare_means()+
#   xlab("")+ylab("Moran index")
# 
# x1 <- moran[setdiff(gene.list[[1]],gene.list[[3]]),]
# x2 <- moran[setdiff(gene.list[[3]],gene.list[[1]]),]
# p2 <- data.frame(name = c(rep("Festem only",nrow(x1)),rep("HVGdisp only",nrow(x2))),
#                  value = c(x1$observed,x2$observed)) %>%
#   ggplot( aes(x=name, y=value, fill=name)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("Festem only" = "red","HVGdisp only" = "#9BA3EB"))+
#   theme_pubr() +
#   theme(
#     legend.position="none",
#     plot.title = element_text(size=11)
#   ) +
#   stat_compare_means()+
#   xlab("")+ylab("Moran index")
# 
# x1 <- moran[setdiff(gene.list[[1]],gene.list[[6]]),]
# x2 <- moran[setdiff(gene.list[[6]],gene.list[[1]]),]
# p3 <- data.frame(name = c(rep("Festem only",nrow(x1)),rep("trendVar only",nrow(x2))),
#            value = c(x1$observed,x2$observed)) %>%
#   ggplot( aes(x=name, y=value, fill=name)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("Festem only" = "red","trendVar only" = "#C0B236"))+
#   theme_pubr() +
#   theme(
#     legend.position="none",
#     plot.title = element_text(size=11)
#   ) +
#   stat_compare_means()+
#   xlab("")+ylab("Moran index")
# 
# pdf("./figures/FigureS12E.pdf",width = 4,height = 12)
# ggarrange(p1,p2,p3,ncol = 1)
# dev.off()

# Figure S14 --------------------------------------------------------------
load("./results/ifnb_ctrl_clustering_tSNE.RData")
ifnb <- readRDS("./results/ifnb_ctrl.rds")
ifnb <- NormalizeData(ifnb)
ifnb <- ifnb[,!ifnb@meta.data$seurat_annotations%in%c("Eryth","Mk")]
marker_list <- c("ISG15","ISG20","IFI6","IFIT1","IFITM3","IFIT3","IFI35",
                 "LY6E","APOBEC3A","RSAD2","OAS1","MX1","TNFSF10","TNFSF13B")
feature_plots <- vector("list",length(marker_list))
for (i in 1:length(feature_plots)){
  data.tmp <- plots.list[[1]]$data
  exp.tmp <- ifnb@assays$RNA@data[marker_list[i],]
  data.tmp <- cbind.data.frame(data.tmp,exp = exp.tmp)
  max.tmp <- max(data.tmp$exp)
  feature_plots[[i]] <-  ggplot(mapping = aes(x=tSNE_1, y=tSNE_2)) + 
    geom_point(data = data.tmp,
               aes(color=exp),
               cex=0.2,shape = 16)+ 
    scale_color_gradient(limits = c(0,max.tmp),low = "grey75",high = "#CD0404")+
    theme_pubr()+theme(legend.position="bottom") +
    labs(title = marker_list[i],color = "IFN score")
}
pdf("./figures/FigureS14.pdf",width = 8,height = 10,onefile = T)
ggarrange(plotlist = feature_plots, nrow = 4, ncol = 4,
          common.legend = T,legend = "bottom")
dev.off()

