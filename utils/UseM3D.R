UseM3D <- function(data, method) {
    log.data <- data/log(2)
    if (method == "M3Drop") {
        st <- system.time({
            #m3drop, returns df with 1 col being qVal. Genes not ranked are set to qVal = 1
            require(M3Drop)
            
            #png(paste0("m3drop.png"))
            norm <- M3DropConvertData(log.data, is.log = TRUE)
            #Res_m3d <- M3DropDifferentialExpression(data,
            #                                        mt_method="fdr", mt_threshold=0.01)
            #dev.off()
            M3Drop_genes <-
                M3DropFeatureSelection(norm,
                                       mt_method = "fdr",
                                       mt_threshold = 1)
            
            re <- data.frame(row.names = rownames(data))
            re$qVal <- 1
            re[match(rownames(M3Drop_genes), rownames(re)), ] <-
                M3Drop_genes$q.value
        })
        
        p.val <- sort(M3Drop_genes$p.value)
        p.adj <- M3Drop_genes$q.value
        names(p.adj) <- as.vector(M3Drop_genes$Gene)
        p.adj <- sort(p.adj)
        genes <- names(p.adj)
        
        
        return(list("p.val" = p.val, "p.adj" = p.adj, "genes" = genes, "time" = st))
    }
    else if (method == "DANB") {
        st <- system.time({
            require(M3Drop)
            require(Matrix)
            counts <- Matrix(log.data, sparse = TRUE)
            count_mat <- NBumiConvertData(counts, is.log = TRUE, is.counts = FALSE)
            DANB_fit <- NBumiFitModel(count_mat)
            
            # Smoothed gene-specific variances
            par(mfrow = c(1, 2))
            stats <- NBumiCheckFitFS(counts = count_mat, fit = DANB_fit, suppress.plot = TRUE)
            print(c(stats$gene_error, stats$cell_error))
            
            NBDropFS <-
                NBumiFeatureSelectionCombinedDrop(
                    DANB_fit,
                    ntop = nrow(log.data),
                    method = "fdr",
                    suppress.plot = TRUE
                )
        })
        
        # res <- list("p.val" = NBDropFS$p.value, "p.adj" = NBDropFS$q.value)
        # names(res$p.val) <- NBDropFS$Gene
        # names(res$p.adj) <- NBDropFS$Gene
        
        p.val <- sort(NBDropFS$p.value)
        p.adj <- NBDropFS$q.value
        names(p.adj) <- as.vector(NBDropFS$Gene)
        p.adj <- sort(p.adj)
        genes <- names(p.adj)
        
        
        return(list("p.val" = p.val, "p.adj" = p.adj, "genes" = genes, "time" = st))
    }
    
}
