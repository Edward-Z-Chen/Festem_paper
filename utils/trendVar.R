#' Returns ordered gene set based on feature selection of highly variable genes using trendVar

trendVarFS <- function(counts, data) {
    st <- system.time({
        sce <- SingleCellExperiment(list(counts = counts, logcounts = data))
        mgvar <- scran::modelGeneVar(x = sce)
        top.hvgs <- scran::getTopHVGs(mgvar, n = nrow(mgvar))
    })
    
    return(list("var.out" = mgvar, "genes" = top.hvgs, "st" = st))
}