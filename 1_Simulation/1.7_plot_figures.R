library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)
calc_fdr <- function(x, nonDE_index, level = 0.05){
  if (sum(x<level,na.rm = T)==0){
    return(0)
  } else {
    sum(x[nonDE_index]<level,na.rm = T)/sum(x<level,na.rm = T)
  }
}
calc_power <- function(x, DE_index, level = 0.05){
  sum(x[DE_index]<0.05,na.rm = T)/length(DE_index)
}

## Colors and shapes for five cell types setting
my_color <- c(rep("#DE2D26",4),rep("#1F78B4",2),rep("#33A02C",2),"#FF7F00",
              "#6A3D9A","#FF99D7",rep("#B15928",2),"grey50")
names(my_color) <- c("Festem_3g","Festem_5g","Festem_7g","Festem_9g",
                     "DESeq2","DESeq2_full","EdgeR","EdgeR_full",
                     "MAST","Wilcoxon","ROSeq",
                     "singleCellHaystack_PCA","singleCellHaystack_UMAP","TN_test")
my_shape <- c(16,15,17,18,16,7,16,7,16,16,16,16,9,16)
names(my_shape) <- names(my_color)
my_color2 <- c(rep("#DE2D26",4),"#6A3D9A","#9BA3EB","#4DAF4A","#FF7F00","#C0B236")
names(my_color2) <- c("Festem_3g","Festem_5g","Festem_7g","Festem_9g","HVGvst",
                      "HVGdisp","DUBStepR","devianceFS","TrendVar")
my_shape2 <- c(16,15,17,18,rep(16,5))
names(my_shape2) <- names(my_color2)


# Figure 2 (A)-(C) & Figure S3 (A.3)-(D.3) --------------------------------
load("./results/NB_400DE_5type_DEG.RData")
adjpvalue.list <- adjpvalue.list[-c(1,10,12,14)]
names(adjpvalue.list)[c(1,9:11)] <- c("Festem_5g",
                                      "Festem_3g",
                                      "Festem_7g",
                                      "Festem_9g")
adjpvalue.list <- adjpvalue.list[c(9,1,10,11,2:8,12,13)]
time.mat <- time.mat[,c(11,2,13,15,3:8,13,16:17)]
peak.memory.usage <- peak.memory.usage[,c(11,2,13,15,3:8,13,16:17)]

tmp <- read.csv("./results/NB_400DE_5celltype_TN_test.csv",header = F,na.strings = "NA")
tmp <- as.matrix(tmp)
tmp <- t(tmp)
tmp <- t(apply(tmp,1,p.adjust,method = "BH"))

adjpvalue.list <- c(adjpvalue.list,list("TN_test" = tmp))
tmp <- read.csv("./results/NB_400DE_5celltype_TN_test_memory.csv",header = F,na.strings = "NA")
peak.memory.usage <- cbind(peak.memory.usage,"TN_test" = tmp[,2])
tmp <- read.csv("./results/NB_400DE_5celltype_TN_test_time.csv",header = F,na.strings = "NA")
time.mat <- cbind(time.mat,"TN_test" = tmp[,1])

colnames(time.mat) <- names(adjpvalue.list)
colnames(peak.memory.usage) <- names(adjpvalue.list)

deg_frame <- data.frame(method = NA,precision = NA,recall = NA,time = NA,index = NA,peak_memory = NA)
for (i in 1:20){
  for (j in 1:length(adjpvalue.list)){
    deg_frame <- rbind.data.frame(deg_frame,
                                  data.frame(method = names(adjpvalue.list)[j],
                                             precision = 1-calc_fdr(adjpvalue.list[[j]][i,],1:19600),
                                             recall = calc_power(adjpvalue.list[[j]][i,],19601:20000),
                                             time = time.mat[i,j]/60,
                                             peak_memory = peak.memory.usage[i,j],
                                             index = i))
  }
}
deg_frame <- deg_frame[-1,]
deg_frame$method <- factor(deg_frame$method,levels = names(adjpvalue.list))
deg_frame <- cbind(deg_frame,
                   F_score = 2/(1/deg_frame$precision + 1/deg_frame$recall))
deg_frame2 <- aggregate(precision ~ method,deg_frame,mean)
deg_frame2 <- cbind(deg_frame2,
                    recall = aggregate(recall ~ method,deg_frame,mean)[,2],
                    F_score = aggregate(F_score ~ method,deg_frame,mean)[,2])

## Figure S3 (C.3) ---------------------------------------------------------
ggplot(data = deg_frame2,mapping = aes(x = recall, y = precision))+
  geom_point(aes(color = method, shape = method),
             size = 4,stroke = 1)+
  # scale_x_continuous(limits = c(0.4,1))+
  scale_y_continuous(limits = c(0,1))+
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
  guides(color = guide_legend(nrow = 4))


## Figure 2 (C) ------------------------------------------------------------
deg_frame2 <- deg_frame2[c(2,5,7,9,10,11,13,14),]
deg_frame2$method <- as.character(deg_frame2$method)
deg_frame2[c(1,7),1] <- c("Festem","singleCellHaystack")
my_color_tmp <- my_color[c(2,5,7,9,10,11,13,14)]
names(my_color_tmp) <- deg_frame2$method

ggplot(data = deg_frame2,mapping = aes(x = recall, y = precision))+
  geom_point(aes(color = method),
             size = 5,stroke = 1)+
  # scale_x_continuous(limits = c(0.4,1))+
  scale_y_continuous(limits = c(0,1))+
  geom_function(fun = ~ {1/(2/0.9-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.8-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.7-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.6-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.5-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.4-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.3-1/.x)},color = "grey75",linetype = "dashed")+
  theme_pubr()+
  geom_hline(yintercept = 1-0.05,linetype = "dashed",color = "red")+
  scale_color_manual(values = my_color_tmp)+
  # scale_shape_manual(values = my_shape)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 2))


## Figure S3 (D.3) ---------------------------------------------------------
time_plot <- apply(time.mat,2,mean,na.rm = T)/60
time_plot <- data.frame(method = names(adjpvalue.list),
                        time = time_plot,
                        memory = apply(peak.memory.usage,2,mean,na.rm = T)/1024)
time_plot$method <- factor(time_plot$method,levels = names(adjpvalue.list))
ggplot(data = time_plot,mapping = aes(x = memory, y = time))+
  geom_point(aes(color = method, shape = method),
             size = 5,stroke = 1)+
  scale_y_log10()+
  # coord_trans(y=double_exp())+
  theme_pubr()+
  scale_color_manual(values = my_color)+
  scale_shape_manual(values = my_shape)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 4))+
  grids("xy",linetype = "dashed",color = "grey75")+
  theme(panel.grid.minor = element_blank())+
  labs(x = "Memory (GiB)", y = "Time (min)")


## Figure S3 (A.3) ---------------------------------------------------------
load("./results/NB_400DE_5type_FS.RData")
FS.ARI.frame <- FS.ARI.frame[,c(1,7:9,2:6)]
colnames(FS.ARI.frame)[1:4] <- paste0("Festem_",c("5g","3g","7g","9g"))
FS.SI.frame <- FS.SI.frame[,c(1,7:9,2:6)]
colnames(FS.SI.frame)[1:4] <- paste0("Festem_",c("5g","3g","7g","9g"))

FS_frame <- apply(FS.ARI.frame,2,mean,na.rm = T)
FS_frame <- data.frame(method = names(FS_frame),
                       ARI = FS_frame,
                       SI = apply(FS.SI.frame,2,mean,na.rm = T))
ggplot(data = FS_frame,mapping = aes(x = SI, y = ARI))+
  geom_point(aes(color = method,shape = method),
             size = 4,stroke = 1)+
  scale_x_continuous(limits = c(0.2,0.85))+
  scale_y_continuous(limits = c(0.5,1))+
  theme_pubr()+
  scale_color_manual(values = my_color2)+
  scale_shape_manual(values = my_shape2)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 2))+
  grids("xy",linetype = "dashed",color = "grey75")+
  theme(panel.grid.minor = element_blank())

## Figure 2 (B) ------------------------------------------------------------  
FS_frame <- FS_frame[c(1,5:9),]
FS_frame$method <- as.character(FS_frame$method)
FS_frame$method[1] <- "Festem"
my_color2_tmp <- my_color2[c(1,5:9)]
names(my_color2_tmp) <- FS_frame$method
ggplot(data = FS_frame,mapping = aes(x = SI, y = ARI))+
  geom_point(aes(color = method),shape = 15,
             size = 5,stroke = 1)+
  scale_x_continuous(limits = c(0.2,0.85))+
  scale_y_continuous(limits = c(0.5,1))+
  theme_pubr()+
  scale_color_manual(values = my_color2_tmp)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 2))+
  grids("xy",linetype = "dashed",color = "grey75")+
  theme(panel.grid.minor = element_blank())

## Figure S3 (B.3) ---------------------------------------------------------
FS.time.mat <- FS.time.mat[,c(1,7:9,2:6)]
colnames(FS.time.mat)[1:4] <- paste0("Festem_",c("5g","3g","7g","9g"))
peak.memory.usage <- peak.memory.usage[,c(2,11,13,15,18:22)]
time_plot <- apply(FS.time.mat,2,mean,na.rm = T)/60
time_plot <- data.frame(method = colnames(FS.time.mat),
                          time = time_plot,
                          memory = apply(peak.memory.usage,2,mean,na.rm = T)/1024)
time_plot$method <- factor(time_plot$method,levels = colnames(FS.time.mat))
ggplot(data = time_plot,mapping = aes(x = memory, y = time))+
  geom_point(aes(color = method, shape = method),
             size = 5,stroke = 1)+
  scale_y_log10()+
  theme_pubr()+
  scale_color_manual(values = my_color2)+
  scale_shape_manual(values = my_shape2)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 4))+
  grids("xy",linetype = "dashed",color = "grey75")+
  theme(panel.grid.minor = element_blank())+
  labs(x = "Memory (GiB)", y = "Time (min)")

## Figure 2 (A) ------------------------------------------------------------
for (i in 1:length(plots.list)){
  data_plot <- plots.list[[i]]$data
  data_plot$cluster <- factor(c(rep(1,1500),rep(2,375),rep(3,375),rep(4,375),rep(5,375)))
  my_color <- RColorBrewer::brewer.pal(8, "Set1")
  names(my_color) <- 1:8
  plots.list[[i]] <- ggplot(data_plot,aes(x=UMAP_1,y=UMAP_2,color=cluster))+
    geom_point(size = 0.01)+theme_void()+
    scale_color_manual(values = my_color)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
    labs(title =  names(FS.label.list)[i])
}
ggarrange(plotlist = plots.list[1:6],nrow = 2,ncol = 3,align = "hv")

# Figure 2 (D)-(F) & Figure S3 (A.4)-(D.4) --------------------------------
load("./results/NB_200DE_5type_DEG.RData")
adjpvalue.list <- adjpvalue.list[-c(1,10,12,14)]
names(adjpvalue.list)[c(1,9:11)] <- c("Festem_5g",
                                      "Festem_3g",
                                      "Festem_7g",
                                      "Festem_9g")
adjpvalue.list <- adjpvalue.list[c(9,1,10,11,2:8,12,13)]
time.mat <- time.mat[,c(11,2,13,15,3:8,13,16:17)]
peak.memory.usage <- peak.memory.usage[,c(11,2,13,15,3:8,13,16:17)]

tmp <- read.csv("./results/NB_400DE_5celltype_TN_test.csv",header = F,na.strings = "NA")
tmp <- as.matrix(tmp)
tmp <- t(tmp)
tmp <- t(apply(tmp,1,p.adjust,method = "BH"))

adjpvalue.list <- c(adjpvalue.list,list("TN_test" = tmp))
tmp <- read.csv("./results/NB_400DE_5celltype_TN_test_memory.csv",header = F,na.strings = "NA")
peak.memory.usage <- cbind(peak.memory.usage,"TN_test" = tmp[,2])
tmp <- read.csv("./results/NB_400DE_5celltype_TN_test_time.csv",header = F,na.strings = "NA")
time.mat <- cbind(time.mat,"TN_test" = tmp[,1])

colnames(time.mat) <- names(adjpvalue.list)
colnames(peak.memory.usage) <- names(adjpvalue.list)

deg_frame <- data.frame(method = NA,precision = NA,recall = NA,time = NA,index = NA,peak_memory = NA)
for (i in 1:20){
  for (j in 1:length(adjpvalue.list)){
    deg_frame <- rbind.data.frame(deg_frame,
                                  data.frame(method = names(adjpvalue.list)[j],
                                             precision = 1-calc_fdr(adjpvalue.list[[j]][i,],1:19800),
                                             recall = calc_power(adjpvalue.list[[j]][i,],19801:20000),
                                             time = time.mat[i,j]/60,
                                             peak_memory = peak.memory.usage[i,j],
                                             index = i))
  }
}
deg_frame <- deg_frame[-1,]
deg_frame$method <- factor(deg_frame$method,levels = names(adjpvalue.list))
deg_frame <- cbind(deg_frame,
                   F_score = 2/(1/deg_frame$precision + 1/deg_frame$recall))
deg_frame2 <- aggregate(precision ~ method,deg_frame,mean)
deg_frame2 <- cbind(deg_frame2,
                    recall = aggregate(recall ~ method,deg_frame,mean)[,2],
                    F_score = aggregate(F_score ~ method,deg_frame,mean)[,2])

## Figure S3 (C.4) ---------------------------------------------------------
ggplot(data = deg_frame2,mapping = aes(x = recall, y = precision))+
  geom_point(aes(color = method, shape = method),
             size = 4,stroke = 1)+
  # scale_x_continuous(limits = c(0.4,1))+
  scale_y_continuous(limits = c(0,1))+
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
  guides(color = guide_legend(nrow = 4))


## Figure 2 (F) ------------------------------------------------------------
deg_frame2 <- deg_frame2[c(2,5,7,9,10,11,13,14),]
deg_frame2$method <- as.character(deg_frame2$method)
deg_frame2[c(1,7),1] <- c("Festem","singleCellHaystack")
my_color_tmp <- my_color[c(2,5,7,9,10,11,13,14)]
names(my_color_tmp) <- deg_frame2$method

ggplot(data = deg_frame2,mapping = aes(x = recall, y = precision))+
  geom_point(aes(color = method),
             size = 5,stroke = 1)+
  # scale_x_continuous(limits = c(0.4,1))+
  scale_y_continuous(limits = c(0,1))+
  geom_function(fun = ~ {1/(2/0.9-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.8-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.7-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.6-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.5-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.4-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.3-1/.x)},color = "grey75",linetype = "dashed")+
  theme_pubr()+
  geom_hline(yintercept = 1-0.05,linetype = "dashed",color = "red")+
  scale_color_manual(values = my_color_tmp)+
  # scale_shape_manual(values = my_shape)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 2))


## Figure S3 (D.4) ---------------------------------------------------------
time_plot <- apply(time.mat,2,mean,na.rm = T)/60
time_plot <- data.frame(method = names(adjpvalue.list),
                        time = time_plot,
                        memory = apply(peak.memory.usage,2,mean,na.rm = T)/1024)
time_plot$method <- factor(time_plot$method,levels = names(adjpvalue.list))
ggplot(data = time_plot,mapping = aes(x = memory, y = time))+
  geom_point(aes(color = method, shape = method),
             size = 5,stroke = 1)+
  scale_y_log10()+
  # coord_trans(y=double_exp())+
  theme_pubr()+
  scale_color_manual(values = my_color)+
  scale_shape_manual(values = my_shape)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 4))+
  grids("xy",linetype = "dashed",color = "grey75")+
  theme(panel.grid.minor = element_blank())+
  labs(x = "Memory (GiB)", y = "Time (min)")


## Figure S3 (A.4) ---------------------------------------------------------
load("./results/NB_200DE_5type_FS.RData")
FS.ARI.frame <- FS.ARI.frame[,c(1,7:9,2:6)]
colnames(FS.ARI.frame)[1:4] <- paste0("Festem_",c("5g","3g","7g","9g"))
FS.SI.frame <- FS.SI.frame[,c(1,7:9,2:6)]
colnames(FS.SI.frame)[1:4] <- paste0("Festem_",c("5g","3g","7g","9g"))

FS_frame <- apply(FS.ARI.frame,2,mean,na.rm = T)
FS_frame <- data.frame(method = names(FS_frame),
                       ARI = FS_frame,
                       SI = apply(FS.SI.frame,2,mean,na.rm = T))
ggplot(data = FS_frame,mapping = aes(x = SI, y = ARI))+
  geom_point(aes(color = method,shape = method),
             size = 5,stroke = 1)+
  # scale_x_continuous(limits = c(0.2,0.85))+
  # scale_y_continuous(limits = c(0.5,1))+
  theme_pubr()+
  scale_color_manual(values = my_color2)+
  scale_shape_manual(values = my_shape2)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 2))+
  grids("xy",linetype = "dashed",color = "grey75")+
  theme(panel.grid.minor = element_blank())

## Figure 2 (E) ------------------------------------------------------------  
FS_frame <- FS_frame[c(1,5:9),]
FS_frame$method <- as.character(FS_frame$method)
FS_frame$method[1] <- "Festem"
my_color2_tmp <- my_color2[c(1,5:9)]
names(my_color2_tmp) <- FS_frame$method
ggplot(data = FS_frame,mapping = aes(x = SI, y = ARI))+
  geom_point(aes(color = method),shape = 15,
             size = 5,stroke = 1)+
  # scale_x_continuous(limits = c(0.2,0.85))+
  # scale_y_continuous(limits = c(0.5,1))+
  theme_pubr()+
  scale_color_manual(values = my_color2_tmp)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 2))+
  grids("xy",linetype = "dashed",color = "grey75")+
  theme(panel.grid.minor = element_blank())

## Figure S3 (B.4) ---------------------------------------------------------
FS.time.mat <- FS.time.mat[,c(1,7:9,2:6)]
colnames(FS.time.mat)[1:4] <- paste0("Festem_",c("5g","3g","7g","9g"))
peak.memory.usage <- peak.memory.usage[,c(2,11,13,15,18:22)]
time_plot <- apply(FS.time.mat,2,mean,na.rm = T)/60
time_plot <- data.frame(method = colnames(FS.time.mat),
                        time = time_plot,
                        memory = apply(peak.memory.usage,2,mean,na.rm = T)/1024)
time_plot$method <- factor(time_plot$method,levels = colnames(FS.time.mat))
ggplot(data = time_plot,mapping = aes(x = memory, y = time))+
  geom_point(aes(color = method, shape = method),
             size = 5,stroke = 1)+
  scale_y_log10()+
  theme_pubr()+
  scale_color_manual(values = my_color2)+
  scale_shape_manual(values = my_shape2)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 4))+
  grids("xy",linetype = "dashed",color = "grey75")+
  theme(panel.grid.minor = element_blank())+
  labs(x = "Memory (GiB)", y = "Time (min)")

## Figure 2 (D) ------------------------------------------------------------
for (i in 1:length(plots.list)){
  data_plot <- plots.list[[i]]$data
  data_plot$cluster <- factor(c(rep(1,1500),rep(2,375),rep(3,375),rep(4,375),rep(5,375)))
  my_color <- RColorBrewer::brewer.pal(8, "Set1")
  names(my_color) <- 1:8
  plots.list[[i]] <- ggplot(data_plot,aes(x=UMAP_1,y=UMAP_2,color=cluster))+
    geom_point(size = 0.01)+theme_void()+
    scale_color_manual(values = my_color)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
    labs(title =  names(FS.label.list)[i])
}
ggarrange(plotlist = plots.list[1:6],nrow = 2,ncol = 3,align = "hv")

## Colors and shapes for two cell types setting
my_color <- c(rep("#DE2D26",4),rep("#1F78B4",2),rep("#33A02C",2),"#FF7F00",
              "#6A3D9A","#FF99D7",rep("#B15928",2),"grey50")
names(my_color) <- c("Festem_2g","Festem_3g","Festem_4g","Festem_5g",
                     "DESeq2","DESeq2_full","EdgeR","EdgeR_full",
                     "MAST","Wilcoxon","ROSeq",
                     "singleCellHaystack_PCA","singleCellHaystack_UMAP","TN_test")
my_shape <- c(16,15,17,18,16,7,16,7,16,16,16,16,9,16)
names(my_shape) <- names(my_color)
my_color2 <- c(rep("#DE2D26",4),"#6A3D9A","#9BA3EB","#4DAF4A","#FF7F00","#C0B236")
names(my_color2) <- c("Festem_2g","Festem_3g","Festem_4g","Festem_5g","HVGvst",
                      "HVGdisp","DUBStepR","devianceFS","TrendVar")
my_shape2 <- c(16,15,17,18,rep(16,5))
names(my_shape2) <- names(my_color2)

# Figure S2 (A)-(C) & Figure S3 (A.1)-(D.1) --------------------------------
load("./results/NB_400DE_2type_DEG.RData")
adjpvalue.list <- adjpvalue.list[-c(1,10,12,14)]
names(adjpvalue.list)[c(1,9:11)] <- c("Festem_2g",
                                      "Festem_3g",
                                      "Festem_4g",
                                      "Festem_5g")
adjpvalue.list <- adjpvalue.list[c(9,1,10,11,2:8,12,13)]
time.mat <- time.mat[,c(11,2,13,15,3:8,13,16:17)]
peak.memory.usage <- peak.memory.usage[,c(11,2,13,15,3:8,13,16:17)]

tmp <- read.csv("./results/NB_400DE_2celltype_TN_test.csv",header = F,na.strings = "NA")
tmp <- as.matrix(tmp)
tmp <- t(tmp)
tmp <- t(apply(tmp,1,p.adjust,method = "BH"))

adjpvalue.list <- c(adjpvalue.list,list("TN_test" = tmp))
tmp <- read.csv("./results/NB_400DE_2celltype_TN_test_memory.csv",header = F,na.strings = "NA")
peak.memory.usage <- cbind(peak.memory.usage,"TN_test" = tmp[,2])
tmp <- read.csv("./results/NB_400DE_2celltype_TN_test_time.csv",header = F,na.strings = "NA")
time.mat <- cbind(time.mat,"TN_test" = tmp[,1])

colnames(time.mat) <- names(adjpvalue.list)
colnames(peak.memory.usage) <- names(adjpvalue.list)

deg_frame <- data.frame(method = NA,precision = NA,recall = NA,time = NA,index = NA,peak_memory = NA)
for (i in 1:20){
  for (j in 1:length(adjpvalue.list)){
    deg_frame <- rbind.data.frame(deg_frame,
                                  data.frame(method = names(adjpvalue.list)[j],
                                             precision = 1-calc_fdr(adjpvalue.list[[j]][i,],1:19600),
                                             recall = calc_power(adjpvalue.list[[j]][i,],19601:20000),
                                             time = time.mat[i,j]/60,
                                             peak_memory = peak.memory.usage[i,j],
                                             index = i))
  }
}
deg_frame <- deg_frame[-1,]
deg_frame$method <- factor(deg_frame$method,levels = names(adjpvalue.list))
deg_frame <- cbind(deg_frame,
                   F_score = 2/(1/deg_frame$precision + 1/deg_frame$recall))
deg_frame2 <- aggregate(precision ~ method,deg_frame,mean)
deg_frame2 <- cbind(deg_frame2,
                    recall = aggregate(recall ~ method,deg_frame,mean)[,2],
                    F_score = aggregate(F_score ~ method,deg_frame,mean)[,2])

## Figure S3 (C.1) ---------------------------------------------------------
ggplot(data = deg_frame2,mapping = aes(x = recall, y = precision))+
  geom_point(aes(color = method, shape = method),
             size = 4,stroke = 1)+
  # scale_x_continuous(limits = c(0.4,1))+
  scale_y_continuous(limits = c(0,1))+
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
  guides(color = guide_legend(nrow = 4))


## Figure S2 (C) ------------------------------------------------------------
deg_frame2 <- deg_frame2[c(2,5,7,9,10,11,13,14),]
deg_frame2$method <- as.character(deg_frame2$method)
deg_frame2[c(1,7),1] <- c("Festem","singleCellHaystack")
my_color_tmp <- my_color[c(2,5,7,9,10,11,13,14)]
names(my_color_tmp) <- deg_frame2$method

ggplot(data = deg_frame2,mapping = aes(x = recall, y = precision))+
  geom_point(aes(color = method),
             size = 5,stroke = 1)+
  # scale_x_continuous(limits = c(0.4,1))+
  scale_y_continuous(limits = c(0,1))+
  geom_function(fun = ~ {1/(2/0.9-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.8-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.7-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.6-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.5-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.4-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.3-1/.x)},color = "grey75",linetype = "dashed")+
  theme_pubr()+
  geom_hline(yintercept = 1-0.05,linetype = "dashed",color = "red")+
  scale_color_manual(values = my_color_tmp)+
  # scale_shape_manual(values = my_shape)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 2))


## Figure S3 (D.1) ---------------------------------------------------------
time_plot <- apply(time.mat,2,mean,na.rm = T)/60
time_plot <- data.frame(method = names(adjpvalue.list),
                        time = time_plot,
                        memory = apply(peak.memory.usage,2,mean,na.rm = T)/1024)
time_plot$method <- factor(time_plot$method,levels = names(adjpvalue.list))
ggplot(data = time_plot,mapping = aes(x = memory, y = time))+
  geom_point(aes(color = method, shape = method),
             size = 5,stroke = 1)+
  scale_y_log10()+
  # coord_trans(y=double_exp())+
  theme_pubr()+
  scale_color_manual(values = my_color)+
  scale_shape_manual(values = my_shape)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 4))+
  grids("xy",linetype = "dashed",color = "grey75")+
  theme(panel.grid.minor = element_blank())+
  labs(x = "Memory (GiB)", y = "Time (min)")


## Figure S3 (A.1) ---------------------------------------------------------
load("./results/NB_400DE_2type_FS.RData")
FS.ARI.frame <- FS.ARI.frame[,c(1,7:9,2:6)]
colnames(FS.ARI.frame)[1:4] <- paste0("Festem_",c("2g","3g","4g","5g"))
FS.SI.frame <- FS.SI.frame[,c(1,7:9,2:6)]
colnames(FS.SI.frame)[1:4] <- paste0("Festem_",c("2g","3g","4g","5g"))

FS_frame <- apply(FS.ARI.frame,2,mean,na.rm = T)
FS_frame <- data.frame(method = names(FS_frame),
                       ARI = FS_frame,
                       SI = apply(FS.SI.frame,2,mean,na.rm = T))
ggplot(data = FS_frame,mapping = aes(x = SI, y = ARI))+
  geom_point(aes(color = method,shape = method),
             size = 5,stroke = 1)+
  # scale_x_continuous(limits = c(0.2,0.85))+
  scale_y_continuous(limits = c(0.7,1))+
  theme_pubr()+
  scale_color_manual(values = my_color2)+
  scale_shape_manual(values = my_shape2)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 2))+
  grids("xy",linetype = "dashed",color = "grey75")+
  theme(panel.grid.minor = element_blank())

## Figure S2 (B) ------------------------------------------------------------  
FS_frame <- FS_frame[c(1,5:9),]
FS_frame$method <- as.character(FS_frame$method)
FS_frame$method[1] <- "Festem"
my_color2_tmp <- my_color2[c(1,5:9)]
names(my_color2_tmp) <- FS_frame$method
ggplot(data = FS_frame,mapping = aes(x = SI, y = ARI))+
  geom_point(aes(color = method),shape = 15,
             size = 5,stroke = 1)+
  # scale_x_continuous(limits = c(0.2,0.85))+
  scale_y_continuous(limits = c(0.7,1))+
  theme_pubr()+
  scale_color_manual(values = my_color2_tmp)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 2))+
  grids("xy",linetype = "dashed",color = "grey75")+
  theme(panel.grid.minor = element_blank())

## Figure S3 (B.1) ---------------------------------------------------------
FS.time.mat <- FS.time.mat[,c(1,7:9,2:6)]
colnames(FS.time.mat)[1:4] <- paste0("Festem_",c("2g","3g","4g","5g"))
peak.memory.usage <- peak.memory.usage[,c(2,11,13,15,18:22)]
time_plot <- apply(FS.time.mat,2,mean,na.rm = T)/60
time_plot <- data.frame(method = colnames(FS.time.mat),
                        time = time_plot,
                        memory = apply(peak.memory.usage,2,mean,na.rm = T)/1024)
time_plot$method <- factor(time_plot$method,levels = colnames(FS.time.mat))
ggplot(data = time_plot,mapping = aes(x = memory, y = time))+
  geom_point(aes(color = method, shape = method),
             size = 5,stroke = 1)+
  scale_y_log10()+
  theme_pubr()+
  scale_color_manual(values = my_color2)+
  scale_shape_manual(values = my_shape2)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 4))+
  grids("xy",linetype = "dashed",color = "grey75")+
  theme(panel.grid.minor = element_blank())+
  labs(x = "Memory (GiB)", y = "Time (min)")

## Figure S2 (A) ------------------------------------------------------------
for (i in 1:length(plots.list)){
  data_plot <- plots.list[[i]]$data
  data_plot$cluster <- factor(c(rep(1,1500),rep(2,375),rep(3,375),rep(4,375),rep(5,375)))
  my_color <- RColorBrewer::brewer.pal(8, "Set1")
  names(my_color) <- 1:8
  plots.list[[i]] <- ggplot(data_plot,aes(x=UMAP_1,y=UMAP_2,color=cluster))+
    geom_point(size = 0.01)+theme_void()+
    scale_color_manual(values = my_color)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
    labs(title =  names(FS.label.list)[i])
}
ggarrange(plotlist = plots.list[1:6],nrow = 2,ncol = 3,align = "hv")

# Figure S2 (D)-(F) & Figure S3 (A.2)-(D.2) --------------------------------
load("./results/NB_200DE_2type_DEG.RData")
adjpvalue.list <- adjpvalue.list[-c(1,10,12,14)]
names(adjpvalue.list)[c(1,9:11)] <- c("Festem_2g",
                                      "Festem_3g",
                                      "Festem_4g",
                                      "Festem_5g")
adjpvalue.list <- adjpvalue.list[c(9,1,10,11,2:8,12,13)]
time.mat <- time.mat[,c(11,2,13,15,3:8,13,16:17)]
peak.memory.usage <- peak.memory.usage[,c(11,2,13,15,3:8,13,16:17)]

tmp <- read.csv("./results/NB_200DE_2celltype_TN_test.csv",header = F,na.strings = "NA")
tmp <- as.matrix(tmp)
tmp <- t(tmp)
tmp <- t(apply(tmp,1,p.adjust,method = "BH"))

adjpvalue.list <- c(adjpvalue.list,list("TN_test" = tmp))
tmp <- read.csv("./results/NB_200DE_2celltype_TN_test_memory.csv",header = F,na.strings = "NA")
peak.memory.usage <- cbind(peak.memory.usage,"TN_test" = tmp[,2])
tmp <- read.csv("./results/NB_200DE_2celltype_TN_test_time.csv",header = F,na.strings = "NA")
time.mat <- cbind(time.mat,"TN_test" = tmp[,1])

colnames(time.mat) <- names(adjpvalue.list)
colnames(peak.memory.usage) <- names(adjpvalue.list)

deg_frame <- data.frame(method = NA,precision = NA,recall = NA,time = NA,index = NA,peak_memory = NA)
for (i in 1:20){
  for (j in 1:length(adjpvalue.list)){
    deg_frame <- rbind.data.frame(deg_frame,
                                  data.frame(method = names(adjpvalue.list)[j],
                                             precision = 1-calc_fdr(adjpvalue.list[[j]][i,],1:19800),
                                             recall = calc_power(adjpvalue.list[[j]][i,],19801:20000),
                                             time = time.mat[i,j]/60,
                                             peak_memory = peak.memory.usage[i,j],
                                             index = i))
  }
}
deg_frame <- deg_frame[-1,]
deg_frame$method <- factor(deg_frame$method,levels = names(adjpvalue.list))
deg_frame <- cbind(deg_frame,
                   F_score = 2/(1/deg_frame$precision + 1/deg_frame$recall))
deg_frame2 <- aggregate(precision ~ method,deg_frame,mean)
deg_frame2 <- cbind(deg_frame2,
                    recall = aggregate(recall ~ method,deg_frame,mean)[,2],
                    F_score = aggregate(F_score ~ method,deg_frame,mean)[,2])

## Figure S3 (C.2) ---------------------------------------------------------
ggplot(data = deg_frame2,mapping = aes(x = recall, y = precision))+
  geom_point(aes(color = method, shape = method),
             size = 5,stroke = 1)+
  # scale_x_continuous(limits = c(0.4,1))+
  scale_y_continuous(limits = c(0,1))+
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
  guides(color = guide_legend(nrow = 4))


## Figure S2 (F) ------------------------------------------------------------
deg_frame2 <- deg_frame2[c(2,5,7,9,10,11,13,14),]
deg_frame2$method <- as.character(deg_frame2$method)
deg_frame2[c(1,7),1] <- c("Festem","singleCellHaystack")
my_color_tmp <- my_color[c(2,5,7,9,10,11,13,14)]
names(my_color_tmp) <- deg_frame2$method

ggplot(data = deg_frame2,mapping = aes(x = recall, y = precision))+
  geom_point(aes(color = method),
             size = 5,stroke = 1)+
  # scale_x_continuous(limits = c(0.4,1))+
  scale_y_continuous(limits = c(0,1))+
  geom_function(fun = ~ {1/(2/0.9-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.8-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.7-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.6-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.5-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.4-1/.x)},color = "grey75",linetype = "dashed")+
  geom_function(fun = ~ {1/(2/0.3-1/.x)},color = "grey75",linetype = "dashed")+
  theme_pubr()+
  geom_hline(yintercept = 1-0.05,linetype = "dashed",color = "red")+
  scale_color_manual(values = my_color_tmp)+
  # scale_shape_manual(values = my_shape)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 2))


## Figure S3 (D.2) ---------------------------------------------------------
time_plot <- apply(time.mat,2,mean,na.rm = T)/60
time_plot <- data.frame(method = names(adjpvalue.list),
                        time = time_plot,
                        memory = apply(peak.memory.usage,2,mean,na.rm = T)/1024)
time_plot$method <- factor(time_plot$method,levels = names(adjpvalue.list))
ggplot(data = time_plot,mapping = aes(x = memory, y = time))+
  geom_point(aes(color = method, shape = method),
             size = 5,stroke = 1)+
  scale_y_log10()+
  # coord_trans(y=double_exp())+
  theme_pubr()+
  scale_color_manual(values = my_color)+
  scale_shape_manual(values = my_shape)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 4))+
  grids("xy",linetype = "dashed",color = "grey75")+
  theme(panel.grid.minor = element_blank())+
  labs(x = "Memory (GiB)", y = "Time (min)")


## Figure S3 (A.2) ---------------------------------------------------------
load("./results/NB_200DE_2type_FS.RData")
FS.ARI.frame <- FS.ARI.frame[,c(1,7:9,2:6)]
colnames(FS.ARI.frame)[1:4] <- paste0("Festem_",c("2g","3g","4g","5g"))
FS.SI.frame <- FS.SI.frame[,c(1,7:9,2:6)]
colnames(FS.SI.frame)[1:4] <- paste0("Festem_",c("2g","3g","4g","5g"))

FS_frame <- apply(FS.ARI.frame,2,mean,na.rm = T)
FS_frame <- data.frame(method = names(FS_frame),
                       ARI = FS_frame,
                       SI = apply(FS.SI.frame,2,mean,na.rm = T))
ggplot(data = FS_frame,mapping = aes(x = SI, y = ARI))+
  geom_point(aes(color = method,shape = method),
             size = 5,stroke = 1)+
  # scale_x_continuous(limits = c(0.2,0.85))+
  # scale_y_continuous(limits = c(0.5,1))+
  theme_pubr()+
  scale_color_manual(values = my_color2)+
  scale_shape_manual(values = my_shape2)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 2))+
  grids("xy",linetype = "dashed",color = "grey75")+
  theme(panel.grid.minor = element_blank())

## Figure S2 (E) ------------------------------------------------------------  
FS_frame <- FS_frame[c(1,5:9),]
FS_frame$method <- as.character(FS_frame$method)
FS_frame$method[1] <- "Festem"
my_color2_tmp <- my_color2[c(1,5:9)]
names(my_color2_tmp) <- FS_frame$method
ggplot(data = FS_frame,mapping = aes(x = SI, y = ARI))+
  geom_point(aes(color = method),shape = 15,
             size = 5,stroke = 1)+
  # scale_x_continuous(limits = c(0.2,0.85))+
  # scale_y_continuous(limits = c(0.7,1))+
  theme_pubr()+
  scale_color_manual(values = my_color2_tmp)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 2))+
  grids("xy",linetype = "dashed",color = "grey75")+
  theme(panel.grid.minor = element_blank())

## Figure S3 (B.2) ---------------------------------------------------------
FS.time.mat <- FS.time.mat[,c(1,7:9,2:6)]
colnames(FS.time.mat)[1:4] <- paste0("Festem_",c("2g","3g","4g","5g"))
peak.memory.usage <- peak.memory.usage[,c(2,11,13,15,18:22)]
time_plot <- apply(FS.time.mat,2,mean,na.rm = T)/60
time_plot <- data.frame(method = colnames(FS.time.mat),
                        time = time_plot,
                        memory = apply(peak.memory.usage,2,mean,na.rm = T)/1024)
time_plot$method <- factor(time_plot$method,levels = colnames(FS.time.mat))
ggplot(data = time_plot,mapping = aes(x = memory, y = time))+
  geom_point(aes(color = method, shape = method),
             size = 5,stroke = 1)+
  scale_y_log10()+
  theme_pubr()+
  scale_color_manual(values = my_color2)+
  scale_shape_manual(values = my_shape2)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0))+
  guides(color = guide_legend(nrow = 4))+
  grids("xy",linetype = "dashed",color = "grey75")+
  theme(panel.grid.minor = element_blank())+
  labs(x = "Memory (GiB)", y = "Time (min)")

## Figure S2 (D) ------------------------------------------------------------
for (i in 1:length(plots.list)){
  data_plot <- plots.list[[i]]$data
  data_plot$cluster <- factor(c(rep(1,1500),rep(2,375),rep(3,375),rep(4,375),rep(5,375)))
  my_color <- RColorBrewer::brewer.pal(8, "Set1")
  names(my_color) <- 1:8
  plots.list[[i]] <- ggplot(data_plot,aes(x=UMAP_1,y=UMAP_2,color=cluster))+
    geom_point(size = 0.01)+theme_void()+
    scale_color_manual(values = my_color)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
    labs(title =  names(FS.label.list)[i])
}
ggarrange(plotlist = plots.list[1:6],nrow = 2,ncol = 3,align = "hv")