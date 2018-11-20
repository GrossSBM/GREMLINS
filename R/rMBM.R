#' Simulates multipartite blocl models
#'

#' @param  v_NQ : number of individual in each Functional Group (FG)
#' @param mat_E : define the architecture of the Multipartite.
#' @param type_int : type of interaction (  adjacency symetric or not or incidence)
#' @param vdistrib : vector of the distributions  (bernoulli or poisson for each network)
#' @param lpi : parameters of the clustering distribution
#' @param ltheta : parameters of the interactions distribution
#' @param namesfg : names of the FG.
#' @return Pseudo-Likelihood, penalty
#' @examples
#' A <- matrix(rbinom(100,1,.2),10,10)
#' type <- "diradj"
#' DefineNetwork(A,"diradj","FG1","FG1")
#' @export




##############################################################################################################
rMBM <- function(v_NQ ,E , type_inter, vdistrib, lpi, ltheta, seed=NULL, namesfg= NULL){


  #####
  n_FG <- length(v_NQ)
  if (is.null(namesfg)) {namesfg <- vapply(1:n_FG,function(i){paste('FG',i,sep = '')},'FG1')}
  if (length(namesfg) != n_FG) {stop("Unmatching number of functional group names")}
  if (length(lpi) != n_FG) {stop("Unmatching size of lpi. Should be of length equal to the number of functional groups")}

  # ord <- order(E[,1])
  # E <- E[ord,]
  # ltheta <- ltheta[ord]
  # vdistrib <- vdistrib[ord]
  # type_inter <- type_inter[ord]


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
