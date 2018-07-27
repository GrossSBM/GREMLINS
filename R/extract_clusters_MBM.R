#' Plot the mesoscopic view  of the estimated MBM
#'
#' #'
#' @param fitted_MBM A fitted Generalized BlockModel
#' @param which.model The index corresponding to the model to plot (default is 1, the best model)
#' @examples
#' npc1 <- 20 # number of nodes per block for functional group 1
#' Q1 <- 3 # blocks   for functional group 1
#' n1 <- npc1 * Q1 # number of nodes for functional group 1
#' Z1 <- diag(Q1)%x%matrix(1,npc1,1)
#' P1 <- matrix(runif(Q1*Q1),Q1,Q1)
#' A <- 1*(matrix(runif(n1*n1),n1,n1)<Z1%*%P1%*%t(Z1)) ## adjacency matrix
#' Agr <- DefineNetwork(A,"diradj","FG1","FG1")  # First network
#' npc2 <- 40 #  number of nodes per block for functional group 2
#' Q2 <- 2 # blocks   for functional group 2
#' n2 <- npc2 * Q2 #  number of nodes for functional group 2
#' Z2 <- diag(Q2)%x%matrix(1,npc2,1)
#' P2 <- matrix(runif(Q1*Q2),Q1,Q2)
#' B <- 1*(matrix(runif(n1*n2),n1,n2)<Z1%*%P2%*%t(Z2)) ## incidence matrix
#' Bgr <- DefineNetwork(B,"inc","FG1","FG2")
#' res <- MultipartiteBM(list(Agr,Bgr),namesFG = NULL,vKmin = 1,vKmax = 10,vKinit = NULL,verbose = TRUE, save=FALSE)
#' extract_clusters_MBM(res,mycol=c('blue','red'))
#' @export

extract_clusters_MBM = function(fitted_MBM,which.model = 1){

  dataR6 <- FormattingData(fitted_MBM$listNet)
  param <- fitted_MBM$fitted.model[[which.model]]$param_estim
  vK_estim <- param$vK
  clusters <- lapply(1:dataR6$Q,function(q){lapply(1:vK_estim[q],function(l){names(param$Z[[q]])[param$Z[[q]]==l]})})
  names(clusters) <- dataR6$namesfg
  return(clusters)
}

