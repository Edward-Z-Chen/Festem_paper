library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(scales)
library(RColorBrewer)
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

# Figure 3 (A) -- top right and Figure S5 second from top ------------------------------------------------
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

# Figure 4 ----------------------------------------------------------------
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

## Figure 4G ----------------------------------------------------------------
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

pdf("./figures/Figure4G.pdf",width = 7,height = 4)

dev.off()
