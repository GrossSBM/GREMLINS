#' Predict NAs in a Collection of Networks from a fitted MBM
#'
#' @param RESMBM a fitted multipartite blockmodel
#' @param whichModel The index corresponding to the model used for prediction (default is 1, the best model)
#' @return the collection of matrices of predictions (probability for binary, intensity for weighted network) a
#' @examples
#' set.seed(302718)
#' n_FG <- 2 #number of functional groups (FG)
#' namesFG <- c('A','B')
#' v_NQ <-  c(60,50) #size of each FG
#' v_K  <- c(3,2) #number of clusters in each functional group
#' list_pi = lapply(1:n_FG,function(q){v = rgamma(v_K[q],1,1); return(v/sum(v))})
#' typeInter <- c( "inc","diradj", "adj")
#' v_distrib <- c('gaussian','bernoulli','poisson')
#' E  <-  rbind(c(1,2),c(2,2),c(1,1))
#' list_theta <- list()
#' list_theta[[1]] <- list()
#' m1 <- rnorm(v_K[E[1,1]] * v_K[E[1,2]],7.5,4 )
#' v1 <- rgamma(v_K[E[1,1]] * v_K[E[1,2]],7.5,4 )
#' list_theta[[1]]$mean  <- matrix(m1,nrow = v_K[E[1,1]], ncol = v_K[E[1,2]] )
#' list_theta[[1]]$var  <-  matrix(v1,nrow = v_K[E[1,1]], ncol = v_K[E[1,2]] )
#' m2 <- rbeta(v_K[E[2,1]] * v_K[E[2,2]],2,2 )
#' list_theta[[2]] <- matrix(m2,nrow = v_K[E[2,1]], ncol = v_K[E[2,2]])
#' m3 <- rgamma(v_K[E[3,1]] * v_K[E[3,2]],6,2 )
#' list_theta[[3]] <- matrix((m3 + t(m3))/2,nrow = v_K[E[3,1]], ncol = v_K[E[3,2]])
#' list_Net <- rMBM(v_NQ ,E , typeInter, v_distrib, list_pi,
#'                 list_theta, namesFG = namesFG)$list_Net
#' res_MBMsimu <- multipartiteBM(list_Net,
#'                               v_distrib,
#'                               namesFG = c('A','B'),
#'                               v_Kinit = c(2,2),nbCores = 2)
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
  theta_mean <- lapply(1:nbNet,function(k){m_k <- switch (v_distrib[k], gaussian = theta[[k]]$mean, theta[[k]])})
  # where are NAs
  Pred <- lapply(1:nbNet, function(k){
        tau[[E[k,1]]]%*%theta_mean[[k]] %*%t(tau[[E[k,2]]])
    })
  return(Pred)

}
