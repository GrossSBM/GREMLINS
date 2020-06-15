#' Predict NAs in a Collection of Networks from a fitted MBM
#'
#' @param RESMBM a fitted multipartite blockmodel
#'
#' @details  works only for Bernoulli distribution on all networks. Be sure to define adjacency matrices with NA on the
#' diagonal if there are no loop.
#' @return the collection of networks where the NAs are replaced with predictions (only for Bernoulli model yet) and
#' when there are NAs an additional field \code{wNA} indicating where are the NAs.
#' @export
#'
#' @examples
predictMBM <- function(RESMBM)
{
  data <- RESMBM$list_Net
  Z <- RESMBM$fittedModel[[1]]$paramEstim$Z
  theta <- RESMBM$fittedModel[[1]]$paramEstim$list_theta

  # where are NAs
  LNetPred <- lapply(1:length(data), function(k){
    L <- data[[k]]
    if (sum(is.na(L$mat))==0) return(L) #no NA then return net
    else {
      wNA <- which(is.na(L$mat),arr.ind = T)
      rowClust <- Z[[L$rowFG]][wNA[,1]]
      colClust <- Z[[L$colFG]][wNA[,2]]
      pred <- L$mat
      pred[wNA] <- theta[[k]][cbind(rowClust,colClust)]
      L$mat <- pred
      L$wNA <- wNA
      return(L)
    }
  })
  return(LNetPred)

}
