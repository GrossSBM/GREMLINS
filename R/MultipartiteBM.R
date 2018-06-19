#' Model selection and estimation of multipartite blockmodels
#'
#' Select the number of blocks per functional group using a stepwise search and estimate parameters
#'
#' @param listNet A list of network (defined via the function DefineNetwork)
#' @param namesFG Names of functional groups (must correspond to names in listNet)
#' @param vKmin A vector of minimal number of blocks per functional group provided in the same order as in namesFG
#' @param vKmax A vector of maximal number of blocks per functional group provided in the same order as in namesFG
#' @param vKinit A vector of initial number of blocks per functional group provided in the same order as in namesFG
#' @param verbose Set to TRUE to display the current step of the search algorithm
#' @param save Set to TRUE to save the estimated parameters for intermediate visited models
#' @return a list of estimated parameters for the different models
#' @examples
#' npc1 <- 20 # nodes per class
#' Q1 <- 3 # classes
#' n1 <- npc1 * Q1 # nodes
#' Z1 <- diag(Q1)%x%matrix(1,npc1,1)
#' P1 <- matrix(runif(Q1*Q1),Q1,Q1)
#' A <- 1*(matrix(runif(n1*n1),n1,n1)<Z1%*%P1%*%t(Z1)) ## adjacency matrix
#' Agr <- DefineNetwork(A,"diradj","FG1","FG1")
#' npc2 <- 30 # nodes per class
#' Q2 <- 2 # classes
#' n2 <- npc2 * Q2 # nodes
#' Z2 <- diag(Q2)%x%matrix(1,npc2,1)
#' P2 <- matrix(runif(Q1*Q2),Q1,Q2)
#' B <- 1*(matrix(runif(n1*n2),n1,n2)<Z1%*%P2%*%t(Z2)) ## incidence matrix
#' Bgr <- DefineNetwork(B,"inc","FG1","FG2")
#' res <- MultipartiteBM(list(Agr,Bgr),namesFG = NULL,vKmin = 1,vKmax = 10,vKinit = NULL,verbose = TRUE, save=FALSE)
#' @export

MultipartiteBM = function(listNet,namesFG = NULL,vKmin = 1,vKmax = 10,vKinit = NULL,verbose = TRUE, save=FALSE)
{



  dataR6 = FormattingData(listNet)
  print("------------Nb of entities in each functional group--------------")
  Nb.entities <- dataR6$v_NQ;
  names(Nb.entities) <- dataR6$namesfg;
  print(Nb.entities)



  #------------------- Check the order of names_FG
  if ( dataR6$Q == 1 ) {namesFG <- dataR6$namefg}

  if ( (length(vKmin) == dataR6$Q) & is.null(namesFG))  {stop("Please specify the names of the Functional Groups")}
  if ( (length(vKmax) == dataR6$Q) & is.null(namesFG))  {stop("Please specify the names of the Functional Groups")}
  if ( (length(vKinit) == dataR6$Q) & is.null(namesFG)) {stop("Please specify the names of the Functional Groups")}


  if ((is.null(namesFG)==FALSE)  & (setequal(namesFG,dataR6$namesfg)==FALSE)) {stop("Unmatching names of Functional Groups")}



  if (is.null(vKmin)) {vKmin <- 1; print("The minimum number of clusters has been set to 1")}
  if (is.null(vKmax)) {vKmax <- 10; print("The maximum number of clusters has been set to 10")}
  if (length(vKmin) == 1) {vKmin <- rep(vKmin,dataR6$Q)}else{if (length(vKmin) != dataR6$Q) {stop("Lower bounds on vK are not of the adequate size")}}
  if (length(vKmax) == 1) {vKmax <- rep(vKmax,dataR6$Q)}else{if (length(vKmax) != dataR6$Q) {stop("Upper bounds on vK are not of the adequate size")}}



  vKmin_permut <- vKmin
  vKmax_permut <- vKmax
  vKinit_permut <- vKinit
  for (q in 1:dataR6$Q) {
    wq <- which(dataR6$namesfg == namesFG[q])
    vKmin_permut[wq] <- vKmin[q]
    vKmax_permut[wq] <- vKmax[q]
    if (length(vKinit) == dataR6$Q) {vKinit_permut[q] <- vKinit[wq]}
  }


  vKinit <- vKinit_permut
  vKmin <- vKmin_permut
  vKmax <- vKmax_permut




  #------------------------  Bounds of the numer of clusters
  if (is.null(vKmin)) {vKmin = 1; print("The minimum number of clusters has been set to 1")}
  if (is.null(vKmax)) {vKmax = 10; print("The maximum number of clusters has been set to 10")}

  if (length(vKmin) == 1) {vKmin = rep(vKmin,dataR6$Q)} else {if (length(vKmin) != dataR6$Q) {stop("Lower bounds on vK are not of the adequate size")}}
  if (length(vKmax) == 1) {vKmax = rep(vKmax,dataR6$Q)} else {if (length(vKmax) != dataR6$Q) {stop("Upper bounds on vK are not of the adequate size")}}

  for (q in 1:dataR6$Q)
  {
    if (vKmax[q]>dataR6$v_NQ[q])
    {
      vKmax[q] = dataR6$v_NQ[q]
      print(paste("Kmax[",q,"]  was set to ",dataR6$v_NQ[q],sep=""))
    }
  }


  if (is.null(vKinit)) {
    vKinit_list <- list(vKmin)
    vKmean <- floor((vKmax + vKmin)/2)
    if (sum(vKmean != vKmin) > 0) { vKinit_list[[2]] <- vKmean }
  }else{vKinit_list <- list(vKinit)}



  ################ ESTIMATION starting from one our two initialisation
  R = dataR6$search_nb_clusters(vKinit_list[[1]],Kmin = vKmin,Kmax = vKmax,verbose = verbose)
  if (length(vKinit_list) > 1) {R <- c(R,dataR6$search_nb_clusters(vKinit_list[[2]],Kmin = vKmin,Kmax = vKmax,verbose = verbose))}


  #-------------------- cleaning the results
  ICL_seq <- sapply(R,function(u){u$ICL})
  o <- order(ICL_seq,decreasing = TRUE)
  R.ordered <- lapply(o,function(i){R[[i]]})
  seq_nb_clust <- cbind(t(sapply(R.ordered,function(u){u$param_estim$vK})),ICL_seq[o],1:length(R))
  seq_nb_clust <- seq_nb_clust[!duplicated(seq_nb_clust[,1:dataR6$Q]),]
  res  <- R.ordered[seq_nb_clust[,dataR6$Q + 2]]


  ############### PRINT RESULT#######################
  if (verbose) {
    mess <- paste(res[[1]]$param_estim$vK,collapse = " " )
    mess <- paste("Best model------ ICL :",round(res[[1]]$ICL,2),". Nb of clusters: ", mess,sep = " ")
    mess <- paste(mess, "for",paste(dataR6$namesfg,collapse = " , " ),"respectively",sep = ' ')
    print(mess)
  }

  ############# RESULTATS #############################"
  if (save) {return(list(list(fitted.model = res ,listNet = listNet)))}else{return(list(fitted.model = list(res[[1]]) ,listNet = listNet))}

}

