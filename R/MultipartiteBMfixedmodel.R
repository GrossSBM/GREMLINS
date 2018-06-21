#' Model selection and estimation of multipartite blockmodels
#'
#' Select the number of blocks per functional group using a stepwise search and estimate parameters
#'
#' @param listNet A list of network (defined via the function DefineNetwork)
#' @param namesFG Names of functional groups (must correspond to names in listNet)
#' @param vK A vector with the numbers of blocks per functional group
#' @param classif.init A list of initial classification for each functional group in the same order as in namesFG
#' @param nb_cores Number of cores used for estimation
#' @return Estimated parameters and a classification
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
#' res <- MultipartiteBMfixedmodel(list(Agr,Bgr),namesFG=c("FG1","FG2"),vK=c(3,2))
#' @export


MultipartiteBMfixedmodel <- function(listNet,namesFG ,vK=NULL,classif.init=NULL,nb_cores=NULL){
  #

  dataR6 = FormattingData(listNet)

  if ( dataR6$Q == 1 ) {namesFG <- dataR6$namefg}

  os <- Sys.info()["sysname"]
  if ((os != 'Windows') & (is.null(nb_cores))) {nb_cores = detectCores(all.tests = FALSE, logical = TRUE) %/% 2}

  Q <- dataR6$Q
  vdistrib <- dataR6$vdistrib


  # Check names FG and permute ----------------------------------------------

  if ((is.null(namesFG)==FALSE)  & (setequal(namesFG,dataR6$namesfg)==FALSE)) {stop("Unmatching names of Functional Groups")}

  if ((is.null(vK))&is.null(classif.init)) {stop("one of vK and classif.init have to be defined")}



  if (!is.null(classif.init))
  {
    vKprov = calc_vK(classif.init)
    if (is.null(vK)) {vK=vKprov}
    else {if (sum(vK!=vKprov)>0) {stop("unconsistent initial classification and vK")}}
  }


  #vK_permut <- vK
  permut_vector = numeric(length(vK))
  for (q in 1:dataR6$Q) {
    permut_vector[q] = which(dataR6$namesfg == namesFG[q])
    #wq <- which(dataR6$namesfg == namesFG[q])
    #vK_permut[q] <- vK[wq]
  }
  vK <- vK[permut_vector]

  if (!is.null(classif.init)) {classif.init = classif.init[permut_vector]}


  #------------------------  Initialisation of the number of clusters and the classification.

  if (is.null(classif.init))
  {
    param.init <- genBMfit$new(vK = vK,vdistrib = vdistrib)
    classif.init <- initialize(dataR6,param.init,method = "CAH")$groups
  }


  #----------------------   Initialisation of the algorithm
  estim.0 <- dataR6$estime(classif.init);
  param.0 <- estim.0$param_estim
  classif.0 <- lapply(1:Q,function(q){max.col(param.0$tau[[q]])})

  classif.0 <- lapply(1:dataR6$Q,function(q){Z_q <- max.col(param.0$tau[[q]]);
  Z_q = match(Z_q, unique(sort(Z_q)))
  names(Z_q) <- dataR6$names_ind[[q]]; return(Z_q)})

  ICL.0 <- estim.0$ICL



  #----------------------------------------------  ALGORITHM



  ### step 1 -> M(+1)
  F_forward_q <-  function(q){
    classif_forward_q <- split_classif(classif.0,q,dataR6,100)
    res <- do.call(c, list(classif_forward_q))
    return(res)}
  list_classif_init_forward = list()
  for (q in 1:dataR6$Q) { list_classif_init_forward <- do.call(c,list(list_classif_init_forward,F_forward_q(q)))}


  if (os=="Windows"){
    all_estim_forward <- lapply(list_classif_init_forward,function(init){
      estim.c.l <- dataR6$estime(init)})
  }else{
    all_estim_forward <- mclapply(list_classif_init_forward,function(init){estim.c.l <- dataR6$estime(init)},mc.cores = nb_cores)
  }

  ### step 2 -> M(-1)
  F_backward_q <-  function(q){
    classif_backward_q <- merge_classif(classif.0,q,1)
    res <- do.call(c, list(classif_backward_q))
    return(res)}
  list_classif_init_backward = list()
  for (q in 1:dataR6$Q) { list_classif_init_backward <- do.call(c,list(list_classif_init_backward,F_backward_q(q)))}


  if (os == "Windows") {
    all_estim_backward <- lapply(list_classif_init_backward,function(init){
    estim.c.l <- dataR6$estime(init)})
  }else{
    all_estim_backward <- mclapply(list_classif_init_backward,function(init){estim.c.l <- dataR6$estime(init)},mc.cores = nb_cores)
  }

  ICL.vec.forward <- vapply(all_estim_forward,function(u){u$ICL},1)
  ICL.vec.backward <- vapply(all_estim_backward,function(u){u$ICL},1)

  w.forward <- which.max(ICL.vec.forward)
  estim.new.forward <- all_estim_forward[[w.forward]]
  param.new.forward <- estim.new.forward$param_estim
  classif.new.forward <- lapply(1:dataR6$Q,function(q){Z_q <- max.col(param.new.forward$tau[[q]]);
    Z_q <- match(Z_q, unique(sort(Z_q)))
    names(Z_q) <- dataR6$names_ind[[q]]; return(Z_q)})



  q_forward <- which(param.new.forward$vK != vK)
  init_forward <- merge_classif(classif.new.forward,q_forward,1)
  if (os == "Windows") {
    last_estim_forward <- lapply(init_forward,function(init){estim.c.l <- dataR6$estime(init)})
  }else{
    last_estim_forward <- mclapply(init_forward,function(init){estim.c.l <- dataR6$estime(init)},mc.cores = nb_cores)
  }


  ###########################""
  w.backward = which.max(ICL.vec.backward)
  estim.new.backward <- all_estim_backward[[w.backward]]
  param.new.backward <- estim.new.backward$param_estim
  classif.new.backward <- lapply(1:dataR6$Q,function(q){Z_q <- max.col(param.new.backward$tau[[q]]);
  Z_q <- match(Z_q, unique(sort(Z_q)))
  names(Z_q) <- dataR6$names_ind[[q]]; return(Z_q)})

  q_backward <- which(param.new.backward$vK != vK)
  init_backward <- split_classif(classif.new.backward,q_backward,dataR6,100)
  if(os=="Windows"){
    last_estim_backward <- lapply(init_backward,function(init){estim.c.l <- dataR6$estime(init)})
  }else{
    last_estim_backward <- mclapply(init_backward,function(init){estim.c.l <- dataR6$estime(init)},mc.cores = nb_cores)
  }


  ICL.vec.forward <- vapply(last_estim_forward,function(u){u$ICL},1)
  ICL.vec.backward <- vapply(last_estim_backward,function(u){u$ICL},1)


  ### now we have find the best of the bests!!!!!
  w.star.1 = which.max(ICL.vec.forward)
  w.star.2 = which.max(ICL.vec.backward)

  best = as.character(which.max(c(ICL.0,ICL.vec.forward[w.star.1],ICL.vec.backward[w.star.2])))

  res <- switch(best,
         "1" = {estim.0},
         "2" = {last_estim_forward[[w.star.1]]},
         "3" = {last_estim_backward[[w.star.2]]})
  res$param_estim$Z <- lapply(1:dataR6$Q,function(q){Z_q <- max.col(res$param_estim$tau[[q]]);
  Z_q <- match(Z_q, unique(sort(Z_q)))
  names(Z_q) <- dataR6$names_ind[[q]]; return(Z_q)})

  res$classif <- lapply(res$param_estim$Z,function(z){sort(z)})
  return(res)
}

###########################"


