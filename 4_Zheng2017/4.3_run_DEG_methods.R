source("./DEG_methods/zheng_deseq2.R")
rm(list = ls())

full_flag <- F
source("./DEG_methods/zheng_edger.R")
rm(list = ls())

full_flag <- T
source("./DEG_methods/zheng_edger.R")
rm(list = ls())

source("./DEG_methods/zheng_mast.R")
rm(list = ls())

source("./DEG_methods/zheng_wilcoxon.R")
rm(list = ls())

source("./DEG_methods/zheng_singleCellHaystack.R")
rm(list = ls())

source("./DEG_methods/zheng_ROSeq.R")
rm(list = ls())