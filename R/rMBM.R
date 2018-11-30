#' Simulate datasets from  multipartite block  models
#'
#' @param  v_NQ : number of individual in each Functional Group (FG)
#' @param  E : define the architecture of the Multipartite.
#' @param typeInter : type of interaction (  adjacency symetric or not or incidence) (vector of size equal to nrow(E) )
#' @param v_distrib : vector of the distributions  (bernoulli or poisson for each network) ( vector of size equal to nrow(E) )
#' @param list_pi  : parameters of the blocks distribution
#' @param list_theta : parameters of the interactions distribution
#' @param seed : set the seed for the random simulation (default value  = NULL)
#' @param namesFG : names of the FG.  (default value  = NULL, then the functional groups are labelled FG1, FG2, etc)
#' @return A list of lists. Each element of the list corresponds to a network : each network is described by a matrix (mat) , a type (diradj, adj, inc), the functional group in row (rowFG) and the functional group in columns (colFG)
#' @examples
#' v_K <- c(3,2,2)
#' n_FG <- 3
#' list_pi <- vector("list", 3);
#' list_pi[[1]] <- c(0.4,0.3,0.3); list_pi[[2]] <- c(0.6,0.4); list_pi[[3]]  <- c(0.6,0.4)
#' E  = rbind(c(1,2),c(2,3),c(2,2))
#' v_distrib <- c('bernoulli','poisson','poisson')
#' typeInter <- c( "inc", "inc"  ,  "adj" )
#' list_theta <- list()
#' list_theta[[1]] <- matrix(rbeta(v_K[E[1,1]] * v_K[E[1,2]],1.5,1.5 ),nrow = v_K[E[1,1]], ncol = v_K[E[1,2]] )
#' list_theta[[2]] <- matrix(rgamma(v_K[E[2,1]] * v_K[E[2,2]],7.5,1 ),nrow = v_K[E[2,1]], ncol = v_K[E[2,2]] )
#' list_theta[[3]] <- matrix(rgamma(v_K[E[3,1]] * v_K[E[3,2]],7.5,1 ),nrow = v_K[E[3,1]], ncol = v_K[E[3,2]] )
#' list_theta[[3]] <- 0.5*(list_theta[[3]] + t(list_theta[[3]])) # symetrisation for network 3
#' v_NQ = c(100,50,40)
#' data_sim <- rMBM(v_NQ ,E , typeInter, v_distrib, list_pi, list_theta, seed=NULL, namesFG= c('A','B','D'))
#' @export

##############################################################################################################
rMBM <- function(v_NQ ,E , typeInter, v_distrib, list_pi, list_theta, seed=NULL, namesFG= NULL, keepClassif = FALSE){

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
  if (checkDistrib != 1) {stop("Unmatching definition of the parameters for the chosen distribution (non in [0,1] if Bernoulli of >0 if Poisson")}


  v_K <- vapply(1:n_FG,function(i){length(list_pi[[i]])},2)
  dataSim <- MBMfit$new(v_K,v_distrib,list_pi = list_pi, list_theta = list_theta)
  dataSim <- dataSim$sim(seed=seed,E,v_NQ,typeInter,keepClassif = keepClassif)

  resSim <- dataSim$networks
  list_net <- lapply(1:n_net,function(i){defineNetwork(resSim$mats[[i]],typeInter  = resSim$typeInter[i],rowFG = namesFG[E[i,1]],colFG = namesFG[E[i,2]])})

  res = list(list_net = list_net)
  if (keepClassif) {res$classif <- dataSim$classif}
  return(res)

}
