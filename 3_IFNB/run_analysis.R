library(getopt)

spec <- matrix(
  c("all_analysis", "a", 0,"logical", "Perform all analysis"
    ),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)


# Download Data -----------------------------------------------------------
library(Seurat)
library(SeuratData)
data_list <- InstalledData()
if (!"ifnb"%in%data_list$Dataset){
  InstallData("ifnb")
}
data("ifnb")
ifnb <- subset(ifnb,subset = orig.ident=="IMMUNE_CTRL")
saveRDS(ifnb,file = "./results/ifnb_ctrl.rds")
data("ifnb")
ifnb <- subset(ifnb,subset = orig.ident!="IMMUNE_CTRL")
saveRDS(ifnb,file = "./results/ifnb_stim.rds")

system("wget --no-check-certificate https://housekeeping.unicamp.br/Housekeeping_GenesHuman.RData")

if (!is.null(opt$all_analysis)){
  system("Rscript ./3.1_preclustering.R")
}
system("Rscript ./3.2_run_Festem.R")
if (!is.null(opt$all_analysis)){
  system("Rscript ./3.3_run_DEG_methods.R")
  system("python 3.4_run_TN_test.py")
}
system("Rscript ./3.5_run_FS_methods.R")
system("Rscript ./3.6_clustering_and_tSNE.R")
system("Rscript ./3.7_calculate_CH_indices.R")
if (!is.null(opt$all_analysis)){
  system("Rscript ./3.8_construct_silver_standard.R")
}

if (!is.null(opt$all_analysis)){
  system("Rscript ./3.10_plot_figures.R -a")
} else{
  system("Rscript ./3.10_plot_figures.R")
}