type_allocate <- function(norm.data, type, sig.level = 0.05, plot_result = F, dispersion = NULL){
  # This function returns a logical vector with the same length 
  # of "type" in which each element denotes whether this gene is 
  # highly expressed in this type
  # norm.data should be normal distributed
  # type is a factor denoting which cell belongs to which type
  # "A" means the highest level
  require(ScottKnott)
  lm.data <- lm(norm.data~0+type)
  SK.result <- SK(lm.data,which = "type",sig.level = sig.level)
  if (plot_result){
    plot(SK.result,dispersion = dispersion)
  }
  clus <- SK.result$out$Result
  clus[,1] <- as.numeric(clus[,1])
  clus <- clus[order(clus[,1],decreasing = T),]
  clus.name <- rownames(clus)
  clus <- clus[,-1]
  if (!is.null(ncol(clus))) clus <- apply(clus, 1, paste0,collapse = "")
  clus <- apply(matrix(clus,nrow = 1), 2, toupper)
  names(clus) <- clus.name
  clus[levels(type)]
}
