#' Simulate datasets from  the multipartite block  model (MBM).
#'
#' \code{rMBM} simulates a collection of networks involving common functional groups of entities. The networks may be directed, undirected or bipartite. The emission distribution of the edges may be Bernoulli, Poisson, Gaussian, Zero-Inflated Gaussian, or Laplace. See the vignette for more information about the model.
#'
#' @param  v_NQ : number of individual in each Functional Group (FG)
#' @param  E : define the architecture of the Multipartite.
#' @param  typeInter : type of interaction in each network: undirected adjacency (adj), directed adjacency (diradj) or incidence (inc).  (vector of size equal to nrow(E) )
#' @param  v_distrib : vector of the distributions: 'bernoulli', 'poisson', 'gaussian', 'ZIgaussian' (for Zero inflated gaussian) or 'laplace'  ( vector of size equal to nrow(E) )
#' @param  list_pi  : parameters of the blocks distribution
#' @param  list_theta : parameters of the interactions distribution. For bernoulli numbers between [0,1], for Poisson positive real number, for Gaussian a list specifying mean and var (plus p0 for ZIgaussian), for Laplace a list with location and scale
#' @param  namesFG : names of the FG.  (default value  = NULL, then the functional groups are labelled FG1, FG2, etc)
#' @param  keepClassif : equal to TRUE if you want to keep the simulated blocks/classification (default value  = FALSE).
#' @param  seed : set the seed for the random simulation (default value  = NULL)
#' @return A list of lists containing the networks (list_net) and if keepClassif = TRUE the classifications (classif)
#'         Each element of  list_net corresponds to a network : each network is a list containing  the matrix (mat) , the type of network(diradj, adj, inc), the functional group in row (rowFG) and the functional group in columns (colFG)
#' @examples
#' namesFG <- c('A','B','C')
#' list_pi = list(c(0.16 ,0.40 ,0.44),c(0.3,0.7),c(0.5,0.5))
#' E  <-  rbind(c(1,2),c(2,3),c(1,1))
#' typeInter <- c( "inc","inc", "adj")
#' v_distrib <- c('ZIgaussian','bernoulli','poisson')
#' list_theta <- list()
#' list_theta[[1]] <- list()
#' list_theta[[1]]$mean  <- matrix(c(6.1, 8.9, 6.6, 9.8, 2.6, 1.0), 3, 2)
#' list_theta[[1]]$var  <-  matrix(c(1.6, 1.6, 1.8, 1.7 ,2.3, 1.5),3, 2)
#' list_theta[[1]]$p0  <-  matrix(c(0.4, 0.1, 0.6, 0.5 , 0.2, 0),3, 2)
#' list_theta[[2]] <- matrix(c(0.7,1.0, 0.4, 0.6),2, 2)
#' m3 <- matrix(c(2.5, 2.6 ,2.2 ,2.2, 2.7 ,3.0 ,3.6, 3.5, 3.3),3,3 )
#' list_theta[[3]] <- (m3 + t(m3))/2
#' dataSim <- rMBM(v_NQ = c(100,50,40) , E = E , typeInter = typeInter,
#'                 v_distrib = v_distrib, list_pi = list_pi,
#'                 list_theta = list_theta, namesFG)
#' list_net <- dataSim$list_Net
#' classifTrue <- dataSim$classif
#' @export

##############################################################################################################
rMBM <- function(v_NQ ,E , typeInter, v_distrib, list_pi, list_theta, namesFG= NULL, keepClassif = FALSE, seed=NULL){

  #####
  n_FG <- length(v_NQ)
  if (is.null(namesFG)) {namesFG <- vapply(1:n_FG,function(i){paste('FG',i,sep = '')},'FG1')}
  if (length(namesFG) != n_FG) {stop("Unmatching number of functional group names")}
  if (length(list_pi) != n_FG) {stop("Unmatching size of list_pi. Should be of length equal to the number of functional groups")}

  ####

  if (length(unique(c(E))) != n_FG) {stop("One or more FG non involved in the networks")}
  if (ncol(E) != 2) {stop("wrong definition of mat_E. mat_E should contain 2 columns")}

  n_net  <- nrow(E)
  if (length(typeInter) != n_net) {stop("Unmatching size of typeInter.  Should be of length equal to the number of rows of mat_E")}
  if (length(v_distrib) != n_net) {stop("Unmatching size of vditstrib. Should be of length equal to the number of rows of mat_E")}
  if (length(list_theta)  !=  n_net) {stop("Unmatching size of list_theta.    Should be of length equal to the number of rows of mat_E")}


  ### verif sur les paramètres des lois d'émission
  checkDistrib  = 1 ;
  if (sum(which(v_distrib == 'bernoulli')) > 0) { checkDistrib <- prod(vapply(which(v_distrib == 'bernoulli'),function(i){ prod(list_theta[[i]] <= 1) *  prod(list_theta[[i]] >= 0) },1))}
  if (sum(v_distrib == 'poisson') > 0) { checkDistrib <- checkDistrib  * prod(vapply(which(v_distrib == 'poisson'),function(i){ prod(list_theta[[i]] > 0) },1))}
  if (checkDistrib != 1) {stop("Unmatching definition of the parameters for the chosen distribution (non in [0,1] if Bernoulli or >0 if Poisson")}


  v_K <- vapply(1:n_FG,function(i){length(list_pi[[i]])},2)
  dataSim <- MBMfit$new(v_K,v_distrib,list_pi = list_pi, list_theta = list_theta)
  dataSim <- dataSim$sim(seed = seed,E,v_NQ,typeInter,keepClassif = keepClassif)

  resSim <- dataSim$networks
  list_Net <- lapply(1:n_net,function(i){defineNetwork(resSim$mats[[i]],typeInter  = resSim$typeInter[i],rowFG = namesFG[E[i,1]],colFG = namesFG[E[i,2]])})

  res = list(list_Net = list_Net)
  if (keepClassif) {res$classif <- dataSim$classif}
  return(res)

}
