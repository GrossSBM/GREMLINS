#' Model selection and estimation of multipartite blockmodels
#'
#' Estimate the parameters and give the clustering for given numbers of blocks
#'
#' @param list_Net A list of network (defined via the function DefineNetwork)
#' @param v_distrib  Type of proababilistic distributions in each network : if 0/1 then bernoulli, if counting then poisson, gaussian or Zero Inflated Gaussian (ZIgaussian) My default  = Bernoulli.
#'                   Must give a vector whose length is the number of networks in list_Net
#' @param namesFG Names of functional groups (must correspond to names in listNet)
#' @param v_K A vector with the numbers of blocks per functional group
#' @param classifInit A list of initial classification for each functional group in the same order as in namesFG
#' @param verbose Set to TRUE to display the current step of the search algorithm
#' @param nbCores Number or cores used for the estimation. Not parallelized on windows. By default : half of the cores
#' @param maxiterVE  Maximum number of iterations in the VE step of the VEM algorithm. Default value  = 1000
#' @param maxiterVEM  Maximum number of iterations of the VEM algorithm. Default value  = 1000

#' @return Estimated parameters and a classification
#' @examples
#' namesFG <- c('A','B')
#' list_pi <- list(c(0.5,0.5),c(0.3,0.7)) # prop of blocks in each FG
#' E  <-  rbind(c(1,2),c(2,2)) # architecture of the multipartite net.
#' typeInter <- c( "inc","diradj")
#' v_distrib <- c('poisson','bernoulli')
#' list_theta <- list()
#' list_theta[[1]]   <- matrix(c(6.1, 8.9, 6.6, 3), 2, 2)
#' list_theta[[2]] <- matrix(c(0.7,1.0, 0.4, 0.6),2, 2)
#' list_Net <- rMBM(v_NQ = c(20,20),E , typeInter, v_distrib, list_pi,
#'                 list_theta, namesFG = namesFG, seed = 2)$list_Net
#' #res_MBMsimu_fixed <- multipartiteBMFixedModel(list_Net, v_distrib,
#' #                                               namesFG = namesFG,
#' #                                               v_K = c(1,2),
#' #                                               nbCores = 2)
#' @export


multipartiteBMFixedModel <- function(list_Net,v_distrib ,namesFG , v_K,  classifInit = NULL, nbCores = NULL, maxiterVE = NULL, maxiterVEM = NULL,verbose = TRUE){




  dataR6 = formattingData(list_Net,v_distrib = v_distrib)
  if ( dataR6$Q == 1 ) {namesFG <- dataR6$namesFG}



  os <- Sys.info()["sysname"]
  if (is.null(nbCores)) {nbCores <- detectCores(all.tests = FALSE, logical = TRUE) %/% 2}

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
  estim0 <- dataR6$estime(classifInit,maxiterVE = maxiterVE, maxiterVEM = maxiterVEM);
  param0 <- estim0$paramEstim
  classif0 <- lapply(1:dataR6$Q,
    function(q){
      Z_q <- max.col(param0$tau[[q]]);
      Z_q = match(Z_q, unique(sort(Z_q)))
      names(Z_q) <- dataR6$namesInd[[q]];
      return(Z_q)}
  )

  #----------------------------------------------  ALGORITHM



  ### step 1 -> M(+1)
  Func_Forward_q <-  function(q){
    classiFunc_Forward_q <- splitClassif(classif0,q,dataR6,100)
    res <- do.call(c, list(classiFunc_Forward_q))
    return(res)}
  list_ClassifInitForward = list()
  for (q in 1:dataR6$Q) { list_ClassifInitForward <- do.call(c,list(list_ClassifInitForward,Func_Forward_q(q)))}





  if (os != 'Windows'){
    if (verbose) {
      mess <- '====================== First Forward Step =================='
      print(mess)
      allEstimForward <- pbmcapply::pbmclapply(list_ClassifInitForward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE , maxiterVEM = maxiterVEM)},mc.cores = nbCores)
    }else{
      allEstimForward <- mclapply(list_ClassifInitForward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE , maxiterVEM = maxiterVEM)},mc.cores = nbCores)
    }
  }else{
    if (verbose) {
      mess <- '====================== First Forward Step =================='
      print(mess)
    }
    L <- length(list_ClassifInitForward)
    cl <- parallel::makeCluster(nbCores)
    parallel::clusterExport(cl, c("dataR6","list_ClassifInitForward", "maxiterVE", "maxiterVEM","L"),envir = environment())
    allEstimForward <- parallel::parLapply(cl, 1:L, function(l){estim.c.l <- dataR6$estime(list_ClassifInitForward[[l]],maxiterVE = maxiterVE, maxiterVEM = maxiterVEM)})
    parallel::stopCluster(cl)
  }



  allEstimForward = dataR6$cleanResults(allEstimForward)
  #----------------------------------- ## step 2 -> M(-1)

  Func_Backward_q <-  function(q){
    classiFunc_Backward_q <- mergeClassif(classif0,q,1)
    res <- do.call(c, list(classiFunc_Backward_q))
    return(res)}
  list_ClassifInitBackward = list()
  for (q in 1:dataR6$Q) { list_ClassifInitBackward <- do.call(c,list(list_ClassifInitBackward,Func_Backward_q(q)))}

  if (os != 'Windows'){
    if (verbose) {
      mess <- '====================== First Backward Step =================='
      print(mess)

      allEstimBackward <- pbmcapply::pbmclapply(list_ClassifInitBackward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE, maxiterVEM = maxiterVEM)},mc.cores = nbCores)
    }else{
      allEstimBackward <- mclapply(list_ClassifInitBackward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE , maxiterVEM = maxiterVEM)},mc.cores = nbCores)
    }
  }else{
    if (verbose) {
      mess <- '====================== First Backward Step =================='
      print(mess)
    }
    L <- length(list_ClassifInitBackward)
    cl <- parallel::makeCluster(nbCores)
    parallel::clusterExport(cl, c("dataR6","list_ClassifInitBackward", "maxiterVE", "maxiterVEM","L"),envir = environment())
    allEstimBackward <- parallel::parLapply(cl, 1:L, function(l){estim.c.l <- dataR6$estime(list_ClassifInitBackward[[l]],maxiterVE = maxiterVE, maxiterVEM = maxiterVEM)})
    parallel::stopCluster(cl)
  }

  allEstimBackward = dataR6$cleanResults(allEstimBackward)


  estimNew.forward <- allEstimForward[[1]]
  paramNew.forward <- estimNew.forward$paramEstim
  classifNew.forward <- lapply(1:dataR6$Q,function(q){Z_q <- max.col(paramNew.forward$tau[[q]]);
  Z_q <- match(Z_q, unique(sort(Z_q)))
  names(Z_q) <- dataR6$namesInd[[q]]; return(Z_q)})



  qForward <- which(paramNew.forward$v_K != v_K)
  initForward <- mergeClassif(classifNew.forward,qForward,1)

  if (os != 'Windows'){
    if (verbose) {
      mess <- '====================== Last Forward Step =================='
      print(mess)
      lastEstimForward <- pbmcapply::pbmclapply(initForward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE , maxiterVEM = maxiterVEM)},mc.cores = nbCores)
    }else{
      lastEstimForward <- mclapply(initForward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE, maxiterVEM = maxiterVEM)},mc.cores = nbCores)
    }
  }else{
    if (verbose) {
      mess <- '====================== Last Forward Step =================='
      print(mess)
    }
    L <- length(initForward)
    cl <- parallel::makeCluster(nbCores)
    parallel::clusterExport(cl, c("dataR6","initForward", "maxiterVE", "maxiterVEM","L"),envir = environment())
    lastEstimForward  <- parallel::parLapply(cl, 1:L, function(l){estim.c.l <- dataR6$estime(initForward[[l]],maxiterVE = maxiterVE, maxiterVEM = maxiterVEM)})
    parallel::stopCluster(cl)
  }

  ###########################""
  estimNewBackward <- allEstimBackward[[1]]
  paramNewBackward <- estimNewBackward$paramEstim
  classifNewBackward <- lapply(1:dataR6$Q,function(q){Z_q <- max.col(paramNewBackward$tau[[q]]);
  Z_q <- match(Z_q, unique(sort(Z_q)))
  names(Z_q) <- dataR6$namesInd[[q]]; return(Z_q)})

  qBackward <- which(paramNewBackward$v_K != v_K)
  initBackward <- splitClassif(classifNewBackward,qBackward,dataR6,100)
  if (os != 'Windows'){
    if (verbose) {
      mess <- '====================== Last Backward Step =================='
      print(mess)
      lastEstimBackward <- pbmcapply::pbmclapply(initBackward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE , maxiterVEM = maxiterVEM)},mc.cores = nbCores)
    }else{
      lastEstimBackward <- mclapply(initBackward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE , maxiterVEM = maxiterVEM)},mc.cores = nbCores)
    }
  }else{
    if (verbose) {
      mess <- '====================== Last Backward Step =================='
      print(mess)
    }
    L <- length(initBackward)
    cl <- parallel::makeCluster(nbCores)
    parallel::clusterExport(cl, c("dataR6","initBackward", "maxiterVE", "maxiterVEM","L"),envir = environment())
    lastEstimBackward <- parallel::parLapply(cl, 1:L, function(l){estim.c.l <- dataR6$estime(initBackward[[l]],maxiterVE = maxiterVE, maxiterVEM = maxiterVEM)})
    parallel::stopCluster(cl)
  }


  allEstim <- vector("list", length =  1 + length(lastEstimForward) + length(lastEstimBackward))
  allEstim[[1]] <- estim0
  for (j in 1:length(lastEstimForward)) {allEstim[[j + 1]] =  lastEstimForward[[j]]}
  for (j in 1:length(lastEstimBackward)) {allEstim[[j + length(lastEstimForward)  + 1  ]] = lastEstimBackward[[j]]}


  allEstim <- dataR6$cleanResults(allEstim)
  res <- allEstim[[1]]

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


