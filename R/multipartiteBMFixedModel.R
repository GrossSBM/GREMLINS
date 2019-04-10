#' Model selection and estimation of multipartite blockmodels
#'
#' Estimate the parameters and give the clustering for given numbers of blocks
#'
#' @param list_Net A list of network (defined via the function DefineNetwork)
#' @param namesFG Names of functional groups (must correspond to names in listNet)
#' @param v_K A vector with the numbers of blocks per functional group
#' @param classifInit A list of initial classification for each functional group in the same order as in namesFG
#' @param nb_cores Number of cores used for estimation
#' @param maxiterVE  Maximum number of iterations if the VE step of the VEM algorithm. By default  = 100
#' @return Estimated parameters and a classification
#' @examples
#' v_K <- c(3,2,2)
#' n_FG <- 3
#' list_pi <- vector("list", 3);
#' list_pi[[1]] <- c(0.4,0.3,0.3); list_pi[[2]] <- c(0.6,0.4); list_pi[[3]]  <- c(0.6,0.4)
#' E  = rbind(c(1,2),c(2,3),c(2,2))
#' v_distrib <- c('bernoulli','poisson','poisson')
#' typeInter <- c( "inc", "inc"  ,  "adj" )
#' list_theta <- list()
#' list_theta[[1]] <- matrix(rbeta(v_K[E[1,1]] * v_K[E[1,2]],1.5,1.5 ),nrow = v_K[E[1,1]], ncol = v_K[E[1,2]])
#' list_theta[[2]] <- matrix(rgamma(v_K[E[2,1]] * v_K[E[2,2]],7.5,1 ),nrow = v_K[E[2,1]], ncol = v_K[E[2,2]])
#' list_theta[[3]] <- matrix(rgamma(v_K[E[3,1]] * v_K[E[3,2]],7.5,1 ),nrow = v_K[E[3,1]], ncol = v_K[E[3,2]])
#' list_theta[[3]] <- 0.5*(list_theta[[3]] + t(list_theta[[3]])) # symetrisation for network 3
#' v_NQ = c(100,50,40)
#' list_Net <- rMBM(v_NQ ,E , typeInter, v_distrib, list_pi, list_theta, seed=NULL, namesFG= c('A','B','D'))$list_Net
#' res <- multipartiteBMFixedModel(list_Net,namesFG = c('A','B','D'),v_distrib = v_distrib,v_K = c(3,2,2))
#' @export


multipartiteBMFixedModel <- function(list_Net,v_distrib ,namesFG , v_K=NULL,  classifInit = NULL, nbCores = NULL, maxiterVE = NULL){




  dataR6 = formattingData(list_Net,v_distrib = v_distrib)
  if ( dataR6$Q == 1 ) {namesFG <- dataR6$nameFG}



  os <- Sys.info()["sysname"]
  if ((os != 'Windows') & (is.null(nbCores))) {nbCores = detectCores(all.tests = FALSE, logical = TRUE) %/% 2}



  # Check names FG and permute ----------------------------------------------
  if ((is.null(namesFG) == FALSE)  & (setequal(namesFG,dataR6$namesFG) == FALSE)) {stop("Unmatching names of Functional Groups")}


  if ((is.null(v_K)) & is.null(classifInit)) {stop("one of v_K or classifInit have to be defined")}
  if (!is.null(classifInit))
  {
    v_Kprov = calcVK(classifInit)
    if (is.null(v_K)) {v_K <-  v_Kprov}
    else {if (sum(v_K != v_Kprov) > 0) {stop("unconsistent initial classification and v_K")}}
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
  estim0 <- dataR6$estime(classifInit,maxiterVE = maxiterVE);
  param0 <- estim0$paramEstim
  classif0 <- lapply(1:dataR6$Q,
    function(q){
      Z_q <- max.col(param0$tau[[q]]);
      Z_q = match(Z_q, unique(sort(Z_q)))
      names(Z_q) <- dataR6$namesInd[[q]];
      return(Z_q)}
  )

  #ICL0 <- estim0$ICL



  #----------------------------------------------  ALGORITHM



  ### step 1 -> M(+1)
  Func_Forward_q <-  function(q){
    classiFunc_Forward_q <- splitClassif(classif0,q,dataR6,100)
    res <- do.call(c, list(classiFunc_Forward_q))
    return(res)}
  list_ClassifInitForward = list()
  for (q in 1:dataR6$Q) { list_ClassifInitForward <- do.call(c,list(list_ClassifInitForward,Func_Forward_q(q)))}


  if (os == "Windows") {
    allEstimForward <- lapply(list_ClassifInitForward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE)})
  }else{
    allEstimForward <- mclapply(list_ClassifInitForward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE)},mc.cores = nbCores)
  }


  allEstimForward = dataR6$cleanResults(allEstimForward)
  #----------------------------------- ## step 2 -> M(-1)

  Func_Backward_q <-  function(q){
    classiFunc_Backward_q <- mergeClassif(classif0,q,1)
    res <- do.call(c, list(classiFunc_Backward_q))
    return(res)}
  list_ClassifInitBackward = list()
  for (q in 1:dataR6$Q) { list_ClassifInitBackward <- do.call(c,list(list_ClassifInitBackward,Func_Backward_q(q)))}


  if (os == "Windows") {
    allEstimBackward <- lapply(list_ClassifInitBackward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE)})
  }else{
    allEstimBackward <- mclapply(list_ClassifInitBackward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE)},mc.cores = nbCores)
  }

  allEstimBackward = dataR6$cleanResults(allEstimBackward)


  estimNew.forward <- allEstimForward[[1]]
  paramNew.forward <- estimNew.forward$paramEstim
  classifNew.forward <- lapply(1:dataR6$Q,function(q){Z_q <- max.col(paramNew.forward$tau[[q]]);
  Z_q <- match(Z_q, unique(sort(Z_q)))
  names(Z_q) <- dataR6$namesInd[[q]]; return(Z_q)})



  qForward <- which(paramNew.forward$v_K != v_K)
  initForward <- mergeClassif(classifNew.forward,qForward,1)
  if (os == "Windows") {
    lastEstimForward <- lapply(initForward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE)})
  }else{
    lastEstimForward <- mclapply(initForward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE)},mc.cores = nbCores)
  }



  ###########################""
  estimNewBackward <- allEstimBackward[[1]]
  paramNewBackward <- estimNewBackward$paramEstim
  classifNewBackward <- lapply(1:dataR6$Q,function(q){Z_q <- max.col(paramNewBackward$tau[[q]]);
  Z_q <- match(Z_q, unique(sort(Z_q)))
  names(Z_q) <- dataR6$namesInd[[q]]; return(Z_q)})

  qBackward <- which(paramNewBackward$v_K != v_K)
  initBackward <- splitClassif(classifNewBackward,qBackward,dataR6,100)
  if (os == "Windows") {
    lastEstimBackward <- lapply(initBackward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE)})
  }else{
    lastEstimBackward <- mclapply(initBackward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE)},mc.cores = nbCores)
  }


  allEstim <- vector("list", length =  1 + length(lastEstimForward) + length(lastEstimBackward))
  allEstim[[1]] <- estim0
  for (j in 1:length(lastEstimForward)) {allEstim[[j + 1]] =  lastEstimForward[[j]]}
  for (j in 1:length(lastEstimBackward)) {allEstim[[j + length(lastEstimForward)  + 1  ]] = lastEstimBackward[[j]]}


  #browser()
  allEstim <- dataR6$cleanResults(allEstim)
  res <- allEstim[[1]]

  #garde t on Z ?
  res$paramEstim$Z <- lapply(1:dataR6$Q,function(q){Z_q <- max.col(res$paramEstim$tau[[q]]);
  Z_q <- match(Z_q, unique(sort(Z_q)))
  names(Z_q) <- dataR6$namesInd[[q]]; return(Z_q)})

  res$classif <- res$paramEstim$Z
  resReturn <- list()
  resReturn$fittedModel <- list(res)
  resReturn$list_Net <- list_Net
  return(resReturn)
}

###########################"


