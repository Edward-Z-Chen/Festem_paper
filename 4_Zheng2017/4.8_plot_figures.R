library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(scales)
library(RColorBrewer)
if (!file.exists("./figures")){
  dir.create("./figures")
}

# Figure 3 (A) -- bottom left and Figure S5 second from bottom ------------------------------------------------
## Summarizing results from various DEG detection methods
pbmc <- readRDS("./results/b1.rds")
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]
pbmc@meta.data$celltype <- factor(pbmc@meta.data$celltype)
results <- matrix(nrow = nrow(pbmc),ncol = 19,
                  dimnames = list(rownames(pbmc),
                                  c("Festem","Festem_gamma0.01","Festem_7g","Festem_13g",
                                    "DEseq2", "DEseq2-full",       
                                    "EdgeR", "EdgeR-full","MAST",              
                                    "Wilcoxon","MAST-f",            
                                    "Wilcoxon-f","FC",
                                    "singleCellHaystack-PCA",
                                    "singleCellHaystack-UMAP",
                                    "singleCellHaystack-PCA-20pc",
                                    "singleCellHaystack-UMAP-20pc","TN_test","ROSeq")))
FC.list <- vector("list",nlevels(pbmc@meta.data$celltype))
for (i in 1:nlevels(pbmc@meta.data$celltype)){
  FC.list[[i]] <- FoldChange(pbmc,ident.1 = levels(pbmc@meta.data$celltype)[i],
                             group.by = "celltype")[,1]
}
FC.list <- matrix(unlist(FC.list),ncol = nlevels(pbmc@meta.data$celltype),byrow = F)
FC <- apply(FC.list,1,function(x){max(abs(x))})
results[,"FC"] <- FC
load("./results/b1_Festem.RData")
results[colnames(em.result),"Festem"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result.9)[em.result.9[2,]<0],"Festem"] <- NA
load("./results/b1_Festem_gamma0.01.RData")
results[colnames(em.result),"Festem_gamma0.01"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result.9)[em.result.9[2,]<0],"Festem_gamma0.01"] <- NA
load("./results/b1_Festem_7g.RData")
results[colnames(em.result),"Festem_7g"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result.9)[em.result.9[2,]<0],"Festem_7g"] <- NA
load("./results/b1_Festem_13g.RData")
results[colnames(em.result),"Festem_13g"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result.9)[em.result.9[2,]<0],"Festem_13g"] <- NA
## DESeq2 failed to run due to memory limit.
load("./results/b1_edgeR.RData")
results[names(EdgeR.result),"EdgeR"] <- p.adjust(EdgeR.result,"BH")
load("./results/b1_edgeR_full.RData")
results[names(EdgeR.result),"EdgeR-full"] <- p.adjust(EdgeR.result,"BH")
load("./results/b1_mast.RData")
results[rownames(mast.results),"MAST"] <- p.adjust(mast.results[,3],"BH")
load("./results/b1_wilcoxon.RData")
wil.result <- matrix(unlist(wil.result),ncol = nlevels(pbmc@meta.data$celltype),byrow = F)
wil.result <- apply(wil.result,1,function(x){min(x)*nlevels(pbmc@meta.data$celltype)})
results[,"Wilcoxon"] <- p.adjust(wil.result,"BH")
results[,"Wilcoxon-f"] <- results[,"Wilcoxon"]
results[results[,"FC"]<=0.2,"Wilcoxon-f"] <- NA
results[,"MAST-f"] <- results[,"MAST"]
results[results[,"FC"]<=0.2,"MAST-f"] <- NA
load("./results/b1_haystack.RData")
results[,"singleCellHaystack-PCA"] <- exp(haystack_pca$results$log.p.adj)
results[,"singleCellHaystack-UMAP"] <- exp(haystack_umap$results$log.p.adj)
load("./results/b1_haystack_20PC.RData")
results[,"singleCellHaystack-PCA-20pc"] <- exp(haystack_pca$results$log.p.adj)
results[,"singleCellHaystack-UMAP-20pc"] <- exp(haystack_umap$results$log.p.adj)
load("./results/b1_ROSeq.RData")
results[,"ROSeq"] <- p.adjust(roseq.tmp,"BH")

tn_test <- read.csv("./results/b1_TN_test.csv",
                    header = F)
tn_test <- as.matrix(tn_test)
tn_test <- apply(tn_test,1,function(x){min(x)*length(x)})
results[,"TN_test"] <- p.adjust(tn_test,"BH")
save(results,file = "./results/b1_DEG_results.RData")

## Figure S5 bottom left -----------------------------------------------------------
calc_house_percent <- function(DEG_list,house_list){
  sum(house_list%in%DEG_list)/length(DEG_list)
}
calc_house_retrieve <- function(DEG_list,house_list){
  sum(DEG_list%in%house_list)/length(house_list)
}
load("./results/b1_silver_standard.RData")
g1 <- rownames(moran_h)[abs(moran_h[,1])<0.005 & moran_h[,2]>0.05]
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
              
names(my_color) <- c("Festem","Festem_gamma0.01","Festem_7g","Festem_13g",
                     "DEseq2","DEseq2-full","EdgeR","EdgeR-full",
                     "MAST","Wilcoxon","MAST-f","Wilcoxon-f","ROSeq",
                     "singleCellHaystack-PCA","singleCellHaystack-UMAP",
                     "singleCellHaystack-PCA-20pc","singleCellHaystack-UMAP-20pc","TN_test")
my_shape <- c(16,15,17,18,16,7,16,7,16,16,9,9,16,10,16,12,13,16)
names(my_shape) <- names(my_color)

pdf("./figures/FigureS5_b1_left.pdf",height = 6,width = 6)
ggplot(data = power_frame2,mapping = aes(x = recall, y = precision))+
  geom_point(aes(color = method, shape = method),
             size = 4,stroke = 1)+
  scale_x_continuous(limits = c(0.4,1))+
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

pdf("./figures/FigureS5_b1_right.pdf",height = 8,width = 8)
ggplot(num_reject,aes(x = method,y = num,fill = method))+
  geom_bar(stat = "identity")+
  theme_pubr()+
  geom_text(aes(label=num), vjust=1.6, color="white", size=3.5)+
  scale_fill_manual(values = my_color)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90))+
  guides(color = guide_legend(nrow = 4))
dev.off()

## Figure 3 (A) -- bottom left ------------------------------------------------
power_frame2 <- power_frame2[c(1,5,7,9,10,14,17,18),]
power_frame2[c(1,6),1] <- c("Festem","singleCellHaystack")
my_color <- my_color[c(1,5,7,9,10,13,15,18)]
names(my_color)[7] <- "singleCellHaystack"
my_shape <- c(16,15,17,18,7,8,9,10)
names(my_shape) <- names(my_color)

pdf("./figures/Figure3A_bottomleft.pdf",height = 4,width = 3.5)
ggplot(data = power_frame2,mapping = aes(x = recall, y = precision))+
  geom_point(aes(color = method, shape = method),
             size = 4,stroke = 1)+
  scale_x_continuous(limits = c(0,1))+
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
num_reject <- filter(num_reject,method!="DEseq2")

pdf("./figures/Figure3A_bottomleft_inset.pdf",height = 4,width = 4)
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

# Figure S7 -- bottom left ---------------------------------------------------
load("./results/b1_hvggenes.RData")
pbmc <- readRDS("./results/b1.rds")
pbmc <- NormalizeData(pbmc)
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]
pbmc <- FindVariableFeatures(pbmc,nfeatures = nrow(pbmc))
hvgvst <- VariableFeatures(pbmc)
pbmc <- FindVariableFeatures(pbmc,nfeatures = nrow(pbmc),selection.method = "disp")
hvgdisp <- VariableFeatures(pbmc)

load("./results/b1_silver_standard.RData")
g1 <- rownames(moran_h)[abs(moran_h[,1])<0.005 & moran_h[,2]>0.05]
moran_nonh <- moran_nonh[!is.na(moran_nonh[,1]),]
g2 <- rownames(moran_nonh)[moran_nonh[,1]>0.1]
g2 <- g2[gene_percent[g2]>0.05]
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
    # require(plyr)
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
pdf("./figures/FigureS7_bottomleft.pdf",width = 6,height = 4)
ggrocs(pROC_list)
dev.off()

# Figure 3 (A) -- bottom right and Figure S5 right ------------------------------------------------
## Summarizing results from various DEG detection methods
pbmc <- readRDS("./results/b2.rds")
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]
pbmc@meta.data$celltype <- factor(pbmc@meta.data$celltype)
results <- matrix(nrow = nrow(pbmc),ncol = 19,
                  dimnames = list(rownames(pbmc),
                                  c("Festem","Festem_gamma0.01","Festem_8g","Festem_12g",
                                    "DEseq2", "DEseq2-full",       
                                    "EdgeR", "EdgeR-full","MAST",              
                                    "Wilcoxon","MAST-f",            
                                    "Wilcoxon-f","FC",
                                    "singleCellHaystack-PCA",
                                    "singleCellHaystack-UMAP",
                                    "singleCellHaystack-PCA-20pc",
                                    "singleCellHaystack-UMAP-20pc","TN_test","ROSeq")))
FC.list <- vector("list",nlevels(pbmc@meta.data$celltype))
for (i in 1:nlevels(pbmc@meta.data$celltype)){
  FC.list[[i]] <- FoldChange(pbmc,ident.1 = levels(pbmc@meta.data$celltype)[i],
                             group.by = "celltype")[,1]
}
FC.list <- matrix(unlist(FC.list),ncol = nlevels(pbmc@meta.data$celltype),byrow = F)
FC <- apply(FC.list,1,function(x){max(abs(x))})
results[,"FC"] <- FC
load("./results/b2_Festem.RData")
results[colnames(em.result),"Festem"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result.9)[em.result.9[2,]<0],"Festem"] <- NA
load("./results/b2_Festem_gamma0.01.RData")
results[colnames(em.result),"Festem_gamma0.01"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result.9)[em.result.9[2,]<0],"Festem_gamma0.01"] <- NA
load("./results/b2_Festem_8g.RData")
results[colnames(em.result),"Festem_8g"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result.9)[em.result.9[2,]<0],"Festem_8g"] <- NA
load("./results/b2_Festem_12g.RData")
results[colnames(em.result),"Festem_12g"] <- p.adjust(em.result[1,],"BH")
results[colnames(em.result.9)[em.result.9[2,]<0],"Festem_12g"] <- NA
## DESeq2 failed to run due to memory limit.
load("./results/b2_edgeR.RData")
results[names(EdgeR.result),"EdgeR"] <- p.adjust(EdgeR.result,"BH")
load("./results/b2_edgeR_full.RData")
results[names(EdgeR.result),"EdgeR-full"] <- p.adjust(EdgeR.result,"BH")
load("./results/b2_mast.RData")
results[rownames(mast.results),"MAST"] <- p.adjust(mast.results[,3],"BH")
load("./results/b2_wilcoxon.RData")
wil.result <- matrix(unlist(wil.result),ncol = nlevels(pbmc@meta.data$celltype),byrow = F)
wil.result <- apply(wil.result,1,function(x){min(x)*nlevels(pbmc@meta.data$celltype)})
results[,"Wilcoxon"] <- p.adjust(wil.result,"BH")
results[,"Wilcoxon-f"] <- results[,"Wilcoxon"]
results[results[,"FC"]<=0.2,"Wilcoxon-f"] <- NA
results[,"MAST-f"] <- results[,"MAST"]
results[results[,"FC"]<=0.2,"MAST-f"] <- NA
load("./results/b2_haystack.RData")
results[,"singleCellHaystack-PCA"] <- exp(haystack_pca$results$log.p.adj)
results[,"singleCellHaystack-UMAP"] <- exp(haystack_umap$results$log.p.adj)
load("./results/b2_haystack_20PC.RData")
results[,"singleCellHaystack-PCA-20pc"] <- exp(haystack_pca$results$log.p.adj)
results[,"singleCellHaystack-UMAP-20pc"] <- exp(haystack_umap$results$log.p.adj)
load("./results/b2_ROSeq.RData")
results[,"ROSeq"] <- p.adjust(roseq.tmp,"BH")

tn_test <- read.csv("./results/b2_TN_test.csv",
                    header = F)
tn_test <- as.matrix(tn_test)
tn_test <- apply(tn_test,1,function(x){min(x)*length(x)})
results[,"TN_test"] <- p.adjust(tn_test,"BH")
save(results,file = "./results/b2_DEG_results.RData")

## Figure S5 bottom left -----------------------------------------------------------
calc_house_percent <- function(DEG_list,house_list){
  sum(house_list%in%DEG_list)/length(DEG_list)
}
calc_house_retrieve <- function(DEG_list,house_list){
  sum(DEG_list%in%house_list)/length(house_list)
}
load("./results/b2_silver_standard.RData")
g1 <- rownames(moran_h)[abs(moran_h[,1])<0.005 & moran_h[,2]>0.05]
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
              
names(my_color) <- c("Festem","Festem_gamma0.01","Festem_8g","Festem_12g",
                     "DEseq2","DEseq2-full","EdgeR","EdgeR-full",
                     "MAST","Wilcoxon","MAST-f","Wilcoxon-f","ROSeq",
                     "singleCellHaystack-PCA","singleCellHaystack-UMAP",
                     "singleCellHaystack-PCA-20pc","singleCellHaystack-UMAP-20pc","TN_test")
my_shape <- c(16,15,17,18,16,7,16,7,16,16,9,9,16,10,16,12,13,16)
names(my_shape) <- names(my_color)

pdf("./figures/FigureS5_b2_left.pdf",height = 6,width = 6)
ggplot(data = power_frame2,mapping = aes(x = recall, y = precision))+
  geom_point(aes(color = method, shape = method),
             size = 4,stroke = 1)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0.5,1))+
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

pdf("./figures/FigureS5_b2_right.pdf",height = 8,width = 8)
ggplot(num_reject,aes(x = method,y = num,fill = method))+
  geom_bar(stat = "identity")+
  theme_pubr()+
  geom_text(aes(label=num), vjust=1.6, color="white", size=3.5)+
  scale_fill_manual(values = my_color)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90))+
  guides(color = guide_legend(nrow = 4))
dev.off()

## Figure 3 (A) -- bottom right ------------------------------------------------
power_frame2 <- power_frame2[c(1,5,7,9,10,14,17,18),]
power_frame2[c(1,6),1] <- c("Festem","singleCellHaystack")
my_color <- my_color[c(1,5,7,9,10,13,15,18)]
names(my_color)[7] <- "singleCellHaystack"
my_shape <- c(16,15,17,18,7,8,9,10)
names(my_shape) <- names(my_color)

pdf("./figures/Figure3A_bottomright.pdf",height = 4,width = 3.5)
ggplot(data = power_frame2,mapping = aes(x = recall, y = precision))+
  geom_point(aes(color = method, shape = method),
             size = 4,stroke = 1)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0.5,1))+
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
num_reject <- filter(num_reject,method!="DEseq2")

pdf("./figures/Figure3A_bottomright_inset.pdf",height = 4,width = 4)
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

# Figure S7 -- bottom right ---------------------------------------------------
load("./results/b2_hvggenes.RData")
pbmc <- readRDS("./results/b2.rds")
pbmc <- NormalizeData(pbmc)
pbmc <- pbmc[,!pbmc@meta.data$celltype%in%c("Hematopoietic stem cell","Megakaryocyte","Plasmacytoid dendritic cell")]
pbmc <- FindVariableFeatures(pbmc,nfeatures = nrow(pbmc))
hvgvst <- VariableFeatures(pbmc)
pbmc <- FindVariableFeatures(pbmc,nfeatures = nrow(pbmc),selection.method = "disp")
hvgdisp <- VariableFeatures(pbmc)

load("./results/b2_silver_standard.RData")
g1 <- rownames(moran_h)[abs(moran_h[,1])<0.005 & moran_h[,2]>0.05]
moran_nonh <- moran_nonh[!is.na(moran_nonh[,1]),]
g2 <- rownames(moran_nonh)[moran_nonh[,1]>0.1]
g2 <- g2[gene_percent[g2]>0.05]
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
    # require(plyr)
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
pdf("./figures/FigureS7_bottomright.pdf",width = 6,height = 4)
ggrocs(pROC_list)
dev.off()
