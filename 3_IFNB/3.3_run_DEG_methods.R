source("./DEG_methods/ifnb_deseq2.R")
rm(list = ls())

full_flag <- F
source("./DEG_methods/ifnb_edger.R")
rm(list = ls())

full_flag <- T
source("./DEG_methods/ifnb_edger.R")
rm(list = ls())

source("./DEG_methods/ifnb_mast.R")
rm(list = ls())

source("./DEG_methods/ifnb_wilcoxon.R")
rm(list = ls())

source("./DEG_methods/ifnb_singleCellHaystack.R")
rm(list = ls())