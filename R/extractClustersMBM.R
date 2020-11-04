#' Extract the clusters in each functional group
#'
#'
#' @param resMBM A fitted Generalized BlockModel
#' @param whichModel The index corresponding to the model to plot (default is 1, the best model)
#' @return a list a length the number of Functional Groups. Each element is a list of length the number of blocks composed of the index of the individuals in each block of each cluster.
#' @export

extractClustersMBM = function(resMBM,whichModel = 1){

  v_distrib <- resMBM$fittedModel[[whichModel]]$paramEstim$v_distrib
  dataR6 <- formattingData(resMBM$list_Net,v_distrib)
  param <- resMBM$fittedModel[[whichModel]]$paramEstim
  vK_estim <- param$v_K
  clusters <- lapply(1:length(vK_estim),function(q){lapply(1:vK_estim[q],function(l){
    namesq <- names(param$Z[[q]])
    if (is.null(namesq)){namesq <- 1:length(param$Z[[q]])}
    clustql <- namesq[param$Z[[q]] == l]
    return(clustql)})})
  names(clusters) <- dataR6$namesFG
  return(clusters)
}

