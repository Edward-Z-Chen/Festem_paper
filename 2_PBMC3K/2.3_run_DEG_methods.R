source("./DEG_methods/pbmc_deseq2.R")
rm(list = ls())

full_flag <- F
source("./DEG_methods/pbmc_edger.R")
rm(list = ls())

full_flag <- T
source("./DEG_methods/pbmc_edger.R")
rm(list = ls())

source("./DEG_methods/pbmc_mast.R")
rm(list = ls())

source("./DEG_methods/pbmc_wilcoxon.R")
rm(list = ls())

source("./DEG_methods/pbmc_singleCellHaystack.R")
rm(list = ls())

source("./DEG_methods/pbmc_ROSeq.R")
rm(list = ls())