#' Model selection and estimation of multipartite blockmodels
#'
#' Estimate the parameters and give the clustering for given numbers of blocks
#'
#' @param list_Net A list of network (defined via the function DefineNetwork)
#' @param v_distrib  Type of proababilistic distributions in each network : if 0/1 then bernoulli, if counting then poisson, gaussian or Zero Inflated Gaussian (ZIgaussian) My default  = Bernoulli.
#'                   Must give a vector whose length is the number of networks in list_Net
#' @param namesFG Names of functional groups (must correspond to names in listNet)
#' @param v_K A vector with the numbers of blocks per functional group
#' @param classifInit A list of initial classifications, one per functional group: either named with the names of the functional groups,
#'               or unnamed and in the same order as in namesFG
#' @param initBM If \code{classifInit = NULL} and \code{initBM = TRUE} (default), the initial classification is obtained by fitting a SBM / LBM
#'               on each network with the package \pkg{blockmodels}, adjusted to \code{v_K} blocks (the best start in terms of ICL is kept).
#'               If \code{initBM = FALSE}, a hierarchical clustering is used. Missing values (NA) are allowed.
#' @param verbose Set to TRUE to display the current step of the search algorithm
#' @param nbCores Number or cores used for the estimation. Not parallelized on windows. By default : half of the cores
#' @param maxiterVE  Maximum number of iterations in the VE step of the VEM algorithm. Default value  = 1000
#' @param maxiterVEM  Maximum number of iterations of the VEM algorithm. Default value  = 1000
#'
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


multipartiteBMFixedModel <- function(list_Net,v_distrib ,namesFG , v_K,  classifInit = NULL, nbCores = NULL, maxiterVE = NULL, maxiterVEM = NULL,verbose = TRUE, initBM = TRUE){




  dataR6 = formattingData(list_Net,v_distrib = v_distrib)
  if ( dataR6$Q == 1 ) {namesFG <- dataR6$namesFG}



  os <- Sys.info()["sysname"]
  if (is.null(nbCores)) {nbCores <-future::availableCores()%/% 2}

  # Check names FG and permute ----------------------------------------------
  if ((is.null(namesFG) == FALSE)  & (setequal(namesFG,dataR6$namesFG) == FALSE)) {stop("Unmatching names of Functional Groups")}
  if ((is.null(v_K)) & is.null(classifInit)) {stop("one of v_K or classifInit have to be defined")}
  if (dataR6$Q > 1 && is.null(namesFG) && (length(v_K) == dataR6$Q)) {stop("Please specify the names of the Functional Groups")}

  # put v_K and classifInit in the internal order of the FGs (dataR6$namesFG)
  if (!is.null(v_K)) {
    if (length(v_K) != dataR6$Q) {stop("v_K must be of length the number of functional groups")}
    v_K <- toInternalOrder(v_K, namesFG, dataR6)
  }
  if (!is.null(classifInit)) {
    classifInit <- formatGivenClassif(classifInit, namesFG, dataR6)
    v_Kprov <- calcVK(classifInit)
    if (is.null(v_K)) {v_K <- v_Kprov}
    else {if (any(v_K != v_Kprov)) {stop("unconsistent initial classification and v_K")}}
  }

  #------------------------  Initialisation of the number of clusters and the classification.
  if (is.null(classifInit)) {
    if (initBM) {
      # blockmodels on each network, adjusted to v_K blocks; the best start (ICL) is kept
      listClassifInit <- initBlockmodels(dataR6, v_Kmin = v_K, v_Kmax = v_K, v_K = v_K, verbose = verbose, nbCores = nbCores)
    } else {
      paramInit <- MBMfit$new(v_K = v_K, v_distrib = v_distrib)
      listClassifInit <- list(initialize(dataR6, paramInit, method = "CAH")$groups)
    }
  } else {
    listClassifInit <- list(classifInit)
  }

  #----------------------   Initialisation of the algorithm
  if (length(listClassifInit) == 1) {
    estim0 <- dataR6$estime(listClassifInit[[1]], maxiterVE = maxiterVE, maxiterVEM = maxiterVEM)
  } else {
    allEstim0 <- future.apply::future_lapply(listClassifInit, function(init) {
      dataR6$estime(init, maxiterVE = maxiterVE, maxiterVEM = maxiterVEM)
    }, future.seed = TRUE)
    allICL0 <- vapply(allEstim0, function(est) if (isTRUE(est$convergence)) est$ICL else -Inf, 1)
    if (all(!is.finite(allICL0))) {allICL0 <- vapply(allEstim0, function(est) est$ICL, 1)}
    estim0 <- allEstim0[[which.max(allICL0)]]
  }
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

  if (verbose) {
      mess <- '====================== First Forward Step =================='
      print(mess)
      progressr::handlers("progress")
      L <- length(list_ClassifInitForward)
      allEstimForward <- progressr::with_progress({
        p <- progressr::progressor(steps = L)
        future.apply::future_lapply(seq_len(L), function(l) {
          p(sprintf("init=%g", l))  # Update progress
          estim.c.l <- dataR6$estime(list_ClassifInitForward[[l]], maxiterVE = maxiterVE, maxiterVEM = maxiterVEM)
          return(estim.c.l)
        })
      })
      }else{
      allEstimForward <- future.apply::future_lapply(list_ClassifInitForward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE , maxiterVEM = maxiterVEM)})
   }
  allEstimForward <- if (length(allEstimForward) > 0) dataR6$cleanResults(allEstimForward) else list()

  #----------------------------------- ## step 2 -> M(-1)

  Func_Backward_q <-  function(q){
    classiFunc_Backward_q <- mergeClassif(classif0,q,1)
    res <- do.call(c, list(classiFunc_Backward_q))
    return(res)}
  list_ClassifInitBackward = list()
  for (q in 1:dataR6$Q) { list_ClassifInitBackward <- do.call(c,list(list_ClassifInitBackward,Func_Backward_q(q)))}


  if (verbose) {
    mess <- '====================== First Backward Step =================='
    print(mess)
    progressr::handlers("progress")
    L <- length(list_ClassifInitBackward)
    allEstimBackward <- progressr::with_progress({
      p <- progressr::progressor(steps = L)
      future.apply::future_lapply(seq_len(L), function(l) {
        p(sprintf("init=%g", l))  # Update progress
        estim.c.l <- dataR6$estime(list_ClassifInitBackward[[l]], maxiterVE = maxiterVE, maxiterVEM = maxiterVEM)
        return(estim.c.l)
      })
    })
  }else{
    allEstimBackward <- future.apply::future_lapply(list_ClassifInitBackward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE , maxiterVEM = maxiterVEM)})
  }
  allEstimBackward <- if (length(allEstimBackward) > 0) dataR6$cleanResults(allEstimBackward) else list()


  # nothing to do if the first step was empty (e.g. K_q = 1 for all q in the backward step)
  lastEstimForward <- list()
  if (length(allEstimForward) > 0) {
    estimNew.forward <- allEstimForward[[1]]
    paramNew.forward <- estimNew.forward$paramEstim
    classifNew.forward <- lapply(1:dataR6$Q,function(q){Z_q <- max.col(paramNew.forward$tau[[q]]);
    Z_q <- match(Z_q, unique(sort(Z_q)))
    names(Z_q) <- dataR6$namesInd[[q]]; return(Z_q)})



    qForward <- which(paramNew.forward$v_K != v_K)
    initForward <- mergeClassif(classifNew.forward,qForward,1)

    if (verbose) {
        mess <- '====================== Last Forward Step =================='
        print(mess)
        progressr::handlers("progress")
        L <- length(initForward)
        lastEstimForward <- progressr::with_progress({
          p <- progressr::progressor(steps = L)
          future.apply::future_lapply(seq_len(L), function(l) {
            p(sprintf("init=%g", l))  # Update progress
            estim.c.l <- dataR6$estime(initForward[[l]], maxiterVE = maxiterVE, maxiterVEM = maxiterVEM)
            return(estim.c.l)
            })
          })
        }else{
        lastEstimForward <-  future.apply::future_lapply(initForward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE, maxiterVEM = maxiterVEM)})
        }
  }

  ###########################""
  # nothing to do if the first step was empty (e.g. K_q = 1 for all q in the backward step)
  lastEstimBackward <- list()
  if (length(allEstimBackward) > 0) {
    estimNewBackward <- allEstimBackward[[1]]
    paramNewBackward <- estimNewBackward$paramEstim
    classifNewBackward <- lapply(1:dataR6$Q,function(q){Z_q <- max.col(paramNewBackward$tau[[q]]);
    Z_q <- match(Z_q, unique(sort(Z_q)))
    names(Z_q) <- dataR6$namesInd[[q]]; return(Z_q)})

    qBackward <- which(paramNewBackward$v_K != v_K)
    initBackward <- splitClassif(classifNewBackward,qBackward,dataR6,100)
    if (verbose) {
        mess <- '====================== Last Backward Step =================='
        print(mess)
        L <- length(initBackward)
        lastEstimBackward <- progressr::with_progress({
          p <- progressr::progressor(steps = L)
          future.apply::future_lapply(seq_len(L), function(l) {
            p(sprintf("init=%g", l))  # Update progress
            estim.c.l <- dataR6$estime(initBackward[[l]], maxiterVE = maxiterVE, maxiterVEM = maxiterVEM)
            return(estim.c.l)
          })
        })
        }else{
        lastEstimBackward <- future.apply::future_lapply(initBackward,function(init){estim.c.l <- dataR6$estime(init, maxiterVE = maxiterVE , maxiterVEM = maxiterVEM)})
    }
  }

  allEstim <- vector("list", length =  1 + length(lastEstimForward) + length(lastEstimBackward))
  allEstim[[1]] <- estim0
  for (j in seq_along(lastEstimForward)) {allEstim[[j + 1]] =  lastEstimForward[[j]]}
  for (j in seq_along(lastEstimBackward)) {allEstim[[j + length(lastEstimForward)  + 1  ]] = lastEstimBackward[[j]]}


  allEstim <- dataR6$cleanResults(allEstim)
  if (length(allEstim) == 0) {
    stop("Convergence was not reached. Increase the number of iterations of the VEM or of the VE")
  }
  res <- allEstim[[1]]

  res$paramEstim$Z <- lapply(1:dataR6$Q,function(q){Z_q <- max.col(res$paramEstim$tau[[q]]);
  Z_q <- match(Z_q, unique(sort(Z_q)))
  names(Z_q) <- dataR6$namesInd[[q]]; return(Z_q)})

  names(res$paramEstim$Z) <- dataR6$namesFG
  names(res$paramEstim$v_K) <- dataR6$namesFG
  res$classif <- res$paramEstim$Z
  resReturn <- list()
  resReturn$fittedModel <- list(res)
  resReturn$list_Net <- list_Net
  return(resReturn)
}

###########################"


