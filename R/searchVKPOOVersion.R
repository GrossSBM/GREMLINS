searchKQ <- function(dataR6, classifInit, Kmin=NULL, Kmax=NULL, nbCores=NULL, verbose = TRUE, maxiterVE = NULL){
  #




  os <- Sys.info()["sysname"]
  if ((os != 'Windows') & (is.null(nbCores))) {nbCores <- detectCores(all.tests = FALSE, logical = TRUE) %/% 2}

  #------
  vKinit = calcVK(classifInit)
  if (length(vKinit) != dataR6$Q) { stop('Length of vKinit incorrect') }
  if (verbose) { print(paste("Searching the numbers of blocks starting from",paste(as.character(vKinit),collapse = " "),"clusters",sep = " "))}


  #----------------------   Initialisation of the algorithm
  niterSearch <- 0;
  ICL.c <- -Inf;
  classifC <- classifInit;
  classifNew <- classifC
  estimNew <- dataR6$estime(classifInit,maxiterVE);


  paramNew <- estimNew$paramEstim
  classifNew <- lapply(1:dataR6$Q,
    function(q){
      Z_q = max.col(paramNew$tau[[q]])
      Z_q = match(Z_q, unique(sort(Z_q)))
      names(Z_q) <- dataR6$namesInd[[q]];
      return(Z_q)
      }
    )

  estimNew$paramEstim$Z = classifNew
  ICLNew <- estimNew$ICL

  if (verbose) {
    mess <- paste(round(c(calcVK(classifNew))),collapse = " " )
    print(paste("ICL :",round(ICLNew,2),". Nb of clusters: ", mess,sep = " "))
  }

  vec.ICL = ICLNew
  RES = list()
  RES[[niterSearch + 1]] <- estimNew;
  #----------------------------------------------  ALGORITHM
  while ((ICLNew > ICL.c) & (niterSearch < 1000)) {
    #
    niterSearch <- niterSearch + 1;

    #
    ICL.c <- ICLNew
    classifC <- classifNew

    # list of new classif deriving from the plitting of one cluster or the merging of 2 clusters in cluterisation classifC
    # (These clusterisations will serve as initisalition of the EM algorithm for models)
    list_classif_init <- sequentialInitialize(classifC ,dataR6,Kmin,Kmax,os);
    L = length(list_classif_init)


    if (os == "Windows") {
      allEstim <- lapply(1:L,function(l){
        #print(l);
        estim.c.l <- dataR6$estime(list_classif_init[[l]],maxiterVE = maxiterVE)})
    }else{
      allEstim <- mclapply(1:L,function(l){estim.c.l <- dataR6$estime(list_classif_init[[l]],maxiterVE = maxiterVE)},mc.cores = nbCores)
    }



    all_estim <- dataR6$cleanResults(allEstim)
    ICLNew <- all_estim[[1]]$ICL
    vec.ICL <- c(vec.ICL,ICLNew)

    if (ICLNew > ICL.c) {

      estimNew <- all_estim[[1]]
      estimNew$paramEstim$Z <- lapply(1:dataR6$Q,function(q){Z_q <- max.col(estimNew$paramEstim$tau[[q]]);
      Z_q = match(Z_q, unique(sort(Z_q)))
      names(Z_q) <- dataR6$namesInd[[q]]; return(Z_q)})

      RES[[niterSearch + 1]] <- estimNew;
      paramNew <- estimNew$paramEstim
      classifNew <- paramNew$Z
      if (verbose) {
        mess <- paste(round(c(calcVK(classifNew))),collapse = " " )
        print(paste("ICL :",round(ICLNew,2),". Nb of clusters: ", mess,sep = " "))
      }
    }

  }

  return(RES)
}

###########################"




