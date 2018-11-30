#' Model selection and estimation of multipartite blockmodels
#'
#' Select the number of blocks per functional group using a stepwise search and estimate parameters
#'
#' @param list_Net A list of network (defined via the function DefineNetwork)
#' @param namesFG Names of functional groups (must correspond to names in listNet)
#' @param vK A vector with the numbers of blocks per functional group
#' @param classifInit A list of initial classification for each functional group in the same order as in namesFG
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
#' Bgr <- defineNetwork(B,"inc","FG1","FG2")
#' res <- multipartiteBMFixedModel(list(Agr,Bgr),namesFG=c("FG1","FG2"),vK=c(3,2))
#' @export


multipartiteBMFixedModel <- function(list_Net,namesFG ,v_K=NULL, classifInit = NULL, nbCores = NULL){
  #

  dataR6 = formattingData(list_Net)


  if ( dataR6$Q == 1 ) {namesFG <- dataR6$nameFG}

  os <- Sys.info()["sysname"]
  if ((os != 'Windows') & (is.null(nbCores))) {nbCores = detectCores(all.tests = FALSE, logical = TRUE) %/% 2}


  v_distrib <- dataR6$v_distrib


  # Check names FG and permute ----------------------------------------------

  if ((is.null(namesFG) == FALSE)  & (setequal(namesFG,dataR6$namesFG) == FALSE)) {stop("Unmatching names of Functional Groups")}

  if ((is.null(v_K)) & is.null(classifInit)) {stop("one of vK and classifInit have to be defined")}



  if (!is.null(classifInit))
  {
    v_Kprov = calcVK(classifInit)
    if (is.null(v_K)) {v_K <-  v_Kprov}
    else {if (sum(v_K != v_Kprov) > 0) {stop("unconsistent initial classification and vK")}}
  }

  permut_vector = numeric(length(v_K))
  for (q in 1:dataR6$Q) {
    permut_vector[q] = which(dataR6$namesFG == namesFG[q])
  }
  v_K <- v_K[permut_vector]

  if (!is.null(classifInit)) {classifInit = classifInit[permut_vector]}


  #------------------------  Initialisation of the number of clusters and the classification.

  if (is.null(classifInit))
  {
    paramInit <- MBMfit$new(v_K = v_K,v_distrib = v_distrib)
    classifInit <- initialize(dataR6,paramInit,method = "CAH")$groups
  }



  #----------------------   Initialisation of the algorithm
  estim0 <- dataR6$estime(classifInit);
  param0 <- estim0$paramEstim
  classif0 <- lapply(1:dataR6$Q,
    function(q){
      Z_q <- max.col(param0$tau[[q]]);
      Z_q = match(Z_q, unique(sort(Z_q)))
      names(Z_q) <- dataR6$namesInd[[q]];
      return(Z_q)}
  )

  ICL0 <- estim0$ICL



  #----------------------------------------------  ALGORITHM



  ### step 1 -> M(+1)
  Func_Forward_q <-  function(q){
    classiFunc_Forward_q <- splitClassif(classif0,q,dataR6,100)
    res <- do.call(c, list(classiFunc_Forward_q))
    return(res)}
  list_ClassifInitForward = list()
  for (q in 1:dataR6$Q) { list_ClassifInitForward <- do.call(c,list(list_ClassifInitForward,Func_Forward_q(q)))}


  if (os == "Windows") {
    allEstimForward <- lapply(list_ClassifInitForward,function(init){estim.c.l <- dataR6$estime(init)})
  }else{
    allEstimForward <- mclapply(list_ClassifInitForward,function(init){estim.c.l <- dataR6$estime(init)},mc.cores = nbCores)
  }


  allEstimForward = dataR6$cleanResults(allEstimForward)
  #----------------------------------- ## step 2 -> M(-1)

  Func_Backward_q <-  function(q){
    classiFunc_Backward_q <- mergeClassif(classif.0,q,1)
    res <- do.call(c, list(classiFunc_Backward_q))
    return(res)}
  list_ClassifInitBackward = list()
  for (q in 1:dataR6$Q) { list_ClassifInitBackward <- do.call(c,list(list_ClassifInitBackward,Func_Backward_q(q)))}


  if (os == "Windows") {
    allEstimBackward <- lapply(list_ClassifInitBackward,function(init){estim.c.l <- dataR6$estime(init)})
  }else{
    allEstimBackward <- mclapply(list_ClassifInitBackward,function(init){estim.c.l <- dataR6$estime(init)},mc.cores = nb_cores)
  }

  allEstimBackward = dataR6$cleanResults(allEstimBackward)


  estimNew.forward <- allEstimForward[[1]]
  paramNew.forward <- estimNew.forward$paramEstim
  classifNew.forward <- lapply(1:dataR6$Q,function(q){Z_q <- max.col(paramNew.forward$tau[[q]]);
  Z_q <- match(Z_q, unique(sort(Z_q)))
  names(Z_q) <- dataR6$namesInd[[q]]; return(Z_q)})



  qForward <- which(paramNew.forward$vK != vK)
  initForward <- mergeClassif(classifNew.forward,qForward,1)
  if (os == "Windows") {
    lastEstimForward <- lapply(initForward,function(init){estim.c.l <- dataR6$estime(init)})
  }else{
    lastEstimForward <- mclapply(initForward,function(init){estim.c.l <- dataR6$estime(init)},mc.cores = nb_cores)
  }



  ###########################""
  estimNewBackward <- allEstimBackward[[1]]
  paramNewBackward <- estimNewBackward$paramEstim
  classifNewBackward <- lapply(1:dataR6$Q,function(q){Z_q <- max.col(paramNewBackward$tau[[q]]);
  Z_q <- match(Z_q, unique(sort(Z_q)))
  names(Z_q) <- dataR6$namesInd[[q]]; return(Z_q)})

  qBackward <- which(paramNewBackward$vK != vK)
  initBackward <- split_classif(classifNewBackward,qBackward,dataR6,100)
  if (os == "Windows") {
    lastEstimBackward <- lapply(initBackward,function(init){estim.c.l <- dataR6$estime(init)})
  }else{
    lastEstimBackward <- mclapply(initBackward,function(init){estim.c.l <- dataR6$estime(init)},mc.cores = nb_cores)
  }


  allEstim <- vector("list", length =  1 + length(lastEstimForward) + length(lastEstimBackward))
  allEstim[[1]] <- estim.0
  for (j in 1:length(lastEstimForward)) {allEstim[[j + 1]] =  lastEstimForward[[j]]}
  for (j in 1:length(lastEstimBackward)) {allEstim[[j + length(lastEstimForward)  + 1  ]] = lastEstimBackward[[j]]}


  #browser()
  allEstim <- dataR6$clean_results(allEstim)
  res <- allEstim[[1]]

  #garde t on Z ?
  res$paramEstim$Z <- lapply(1:dataR6$Q,function(q){Z_q <- max.col(res$paramEstim$tau[[q]]);
  Z_q <- match(Z_q, unique(sort(Z_q)))
  names(Z_q) <- dataR6$namesInd[[q]]; return(Z_q)})

  res$classif <- res$paramEstim$Z
  return(res)
}

###########################"


