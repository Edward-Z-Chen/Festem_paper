source("./DEG_methods/b1_deseq2.R")
rm(list = ls())

full_flag <- F
source("./DEG_methods/b1_edger.R")
rm(list = ls())

full_flag <- T
source("./DEG_methods/b1_edger.R")
rm(list = ls())

source("./DEG_methods/b1_mast.R")
rm(list = ls())

source("./DEG_methods/b1_wilcoxon.R")
rm(list = ls())

source("./DEG_methods/b1_singleCellHaystack.R")
rm(list = ls())

source("./DEG_methods/b1_ROSeq.R")
rm(list = ls())