#' Predict NAs in a Collection of Networks from a fitted MBM
#'
#' @param RESMBM a fitted multipartite blockmodel
#' @param whichModel The index corresponding to the model used for prediction (default is 1, the best model)
#' @return the collection of matrices of predictions (probability for binary, intensity for weighted network) a
#' @examples
#' namesFG <- c('A','B')
#' list_pi <- list(c(0.5,0.5),c(0.3,0.7)) # prop of blocks in each FG
#' E  <-  rbind(c(1,2),c(2,2)) # architecture of the multipartite net.
#' typeInter <- c( "inc","diradj")
#' v_distrib <- c('gaussian','bernoulli')
#' list_theta <- list()
#' list_theta[[1]] <- list()
#' list_theta[[1]]$mean  <- matrix(c(6.1, 8.9, 6.6, 3), 2, 2)
#' list_theta[[1]]$var  <-  matrix(c(1.6, 1.6, 1.8, 1.5),2, 2)
#' list_theta[[2]] <- matrix(c(0.7,1.0, 0.4, 0.6),2, 2)
#' list_Net <- rMBM(v_NQ = c(30,30),E , typeInter, v_distrib, list_pi,
#'                 list_theta, namesFG = namesFG, seed = 2)$list_Net
#' res_MBMsimu <- multipartiteBM(list_Net, v_distrib,
#'                               namesFG = c('A','B'), v_Kinit = c(2,2),
#'                               nbCores = 2,initBM = FALSE)
#' pred <- predictMBM(res_MBMsimu)
#' @export
#'
predictMBM <- function(RESMBM,whichModel = 1)
{
  nbNet <- length(RESMBM$list_Net)
  res <- RESMBM$fittedModel[[whichModel]]$paramEstim
  v_distrib <- res$v_distrib
  E <- formattingData(RESMBM$list_Net,v_distrib)$E
  tau <- res$tau

  theta <-  res$list_theta
  theta_mean <- lapply(1:nbNet,function(k){m_k <- switch (v_distrib[k],
                                                          gaussian = theta[[k]]$mean,
                                                          ZIgaussian = ((1-theta[[k]]$p0)>0.5)*theta[[k]]$mean,
                                                          theta[[k]])})
  # where are NAs
  Pred <- lapply(1:nbNet, function(k){
        tau[[E[k,1]]]%*%theta_mean[[k]] %*%t(tau[[E[k,2]]])
    })
  return(Pred)

}
