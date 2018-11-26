#' Simulate datasets from  multipartite block  models
#'
#' @param  v_NQ : number of individual in each Functional Group (FG)
#' @param  E : define the architecture of the Multipartite.
#' @param type_inter : type of interaction (  adjacency symetric or not or incidence) (vector of size equal to nrow(E) )
#' @param vdistrib : vector of the distributions  (bernoulli or poisson for each network) ( vector of size equal to nrow(E) )
#' @param lpi  : parameters of the blocks distribution
#' @param ltheta : parameters of the interactions distribution
#' @param seed : set the seed for the random simulation (default value  = NULL)
#' @param namesfg : names of the FG.  (default value  = NULL, then the functional groups are labelled FG1, FG2, etc)
#' @return A list of lists. Each element of the list corresponds to a network : each network is described by a matrix (mat) , a type (diradj, adj, inc), the functional group in row (rowFG) and the functional group in columns (colFG)
#' @examples
#' v_K <- c(3,2,2)
#' n_FG <- 3
#' lpi <- vector("list", 3);
#' lpi[[1]] <- c(0.4,0.3,0.3); lpi[[2]] <- c(0.6,0.4); lpi[[3]]  <- c(0.6,0.4)
#' E  = rbind(c(1,2),c(2,3),c(2,2))
#' vdistrib <- c('bernoulli','poisson','poisson')
#' type_inter <- c( "inc", "inc"  ,  "adj" )
#' ltheta <- list()
#' ltheta[[1]] <- matrix(rbeta(v_K[E[1,1]] * v_K[E[1,2]],1.5,1.5 ),nrow = v_K[E[1,1]], ncol = v_K[E[1,2]] )
#' ltheta[[2]] <- matrix(rgamma(v_K[E[2,1]] * v_K[E[2,2]],7.5,1 ),nrow = v_K[E[2,1]], ncol = v_K[E[2,2]] )
#' ltheta[[3]] <- matrix(rgamma(v_K[E[3,1]] * v_K[E[3,2]],7.5,1 ),nrow = v_K[E[3,1]], ncol = v_K[E[3,2]] )
#' ltheta[[3]] <- 0.5*(ltheta[[3]] + t(ltheta[[3]])) # symetrisation for network 3
#' v_NQ = c(100,50,40)
#' data_sim <- rMBM(v_NQ ,E , type_inter, vdistrib, lpi, ltheta, seed=NULL, namesfg= c('A','B','D'))
#' @export

##############################################################################################################
rMBM <- function(v_NQ ,E , type_inter, vdistrib, lpi, ltheta, seed=NULL, namesfg= NULL){

  #####
  n_FG <- length(v_NQ)
  if (is.null(namesfg)) {namesfg <- vapply(1:n_FG,function(i){paste('FG',i,sep = '')},'FG1')}
  if (length(namesfg) != n_FG) {stop("Unmatching number of functional group names")}
  if (length(lpi) != n_FG) {stop("Unmatching size of lpi. Should be of length equal to the number of functional groups")}

  ####

  if (length(unique(c(E))) != n_FG) {stop("One or more FG non involved in the networks")}
  if (ncol(E) != 2) {stop("wrong definition of mat_E. mat_E should contain 2 columns")}

  n_net  <- nrow(E)
  if (length(type_inter) != n_net) {stop("Unmatching size of type_inter.  Should be of length equal to the number of rows of mat_E")}
  if (length(vdistrib) != n_net) {stop("Unmatching size of vditstrib. Should be of length equal to the number of rows of mat_E")}
  if (length(ltheta)  !=  n_net) {stop("Unmatching size of ltheta.    Should be of length equal to the number of rows of mat_E")}




  ### verif sur les paramètres des lois d'émission
  check_distrib  = 1 ;
  if (sum(which(vdistrib == 'bernoulli')) > 0) { check_distrib <- prod(vapply(which(vdistrib == 'bernoulli'),function(i){ prod(ltheta[[i]] <= 1) *  prod(ltheta[[i]] >= 0) },1))}
  if (sum(vdistrib == 'poisson') > 0) { check_distrib <- check_distrib  * prod(vapply(which(vdistrib == 'poisson'),function(i){ prod(ltheta[[i]] > 0) },1))}
  if (check_distrib != 1) {stop("Unmatching definition of the parameters for the chosen distribution (non in [0,1] if Bernoulli of >0 if Poisson")}


  vK <- vapply(1:n_FG,function(i){length(lpi[[i]])},2)
  data_sim <- genBMfit$new(vK,vdistrib,lpi = lpi, ltheta = ltheta)
  res_sim <- data_sim$sim(seed=seed,E,v_NQ,type_inter)
  list_net <- lapply(1:n_net,function(i){DefineNetwork(res_sim$mats[[i]],type  = res_sim$type_inter[i],rowFG = namesfg[E[i,1]],colFG = namesfg[E[i,2]])})

  return(list_net)

}
