#' Predict NAs in a Collection of Networks from a fitted MBM
#'
#' @param RESMBM a fitted multipartite blockmodel
#' @param whichModel The index corresponding to the model used for prediction (default is 1, the best model)
#' @return the collection of matrices of predictions (probability for binary, intensity for weighted network) a
#' @examples
#' namesFG <- c('A','B')
#' list_pi = list(c(0.16 ,0.40 ,0.44),c(0.3,0.7))
#' E  <-  rbind(c(1,2),c(2,2),c(1,1))
#' typeInter <- c( "inc","diradj", "adj")
#' v_distrib <- c('gaussian','bernoulli','poisson')
#' list_theta <- list()
#' list_theta[[1]] <- list()
#' list_theta[[1]]$mean  <- matrix(c(6.1, 8.9, 6.6, 9.8, 2.6, 1.0), 3, 2)
#' list_theta[[1]]$var  <-  matrix(c(1.6, 1.6, 1.8, 1.7 ,2.3, 1.5),3, 2)
#' list_theta[[2]] <- matrix(c(0.7,1.0, 0.4, 0.6),2, 2)
#' m3 <- matrix(c(2.5, 2.6 ,2.2 ,2.2, 2.7 ,3.0 ,3.6, 3.5, 3.3),3,3 )
#' list_theta[[3]] <- (m3 + t(m3))/2
#' list_Net <- rMBM(v_NQ = c(60,50),E , typeInter, v_distrib, list_pi,
#'                 list_theta, namesFG = namesFG, seed = 2)$list_Net
#'res_MBMsimu <- multipartiteBM(list_Net, v_distrib,
#'                              namesFG = c('A','B'), v_Kinit = c(2,2),
#'                              nbCores = 2,initBM = FALSE)
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
                                                          ZIgaussian = (1-theta[[k]]$p0)*theta[[k]]$mean,
                                                          theta[[k]])})
  # where are NAs
  Pred <- lapply(1:nbNet, function(k){
        tau[[E[k,1]]]%*%theta_mean[[k]] %*%t(tau[[E[k,2]]])
    })
  return(Pred)

}
