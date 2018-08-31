search_KQ <- function(data,classif.init,Kmin=NULL,Kmax=NULL,nb_cores=NULL,verbose = TRUE){
  #



  os <- Sys.info()["sysname"]
  if ((os != 'Windows') & (is.null(nb_cores))) {nb_cores <- detectCores(all.tests = FALSE, logical = TRUE) %/% 2}
  vdistrib <- data$vdistrib


  #------


  # if (length(vKinit) == 1)
  # {
  #   vKinit_temp <- switch(vKinit,
  #          min =  {Kmin},
  #          mean = {floor((Kmax + Kmin)/2)})
  #   if  (is.numeric(vKinit))
  #     vKinit=rep(vKinit,data$Q)
  #   else vKinit = vKinit_temp
  # }

  vKinit = calc_vK(classif.init)

    if (length(vKinit) != data$Q){ stop('Length of vKinit incorrect') }

    if (verbose) { print(paste("Searching the numbers of blocks starting from",paste(as.character(vKinit),collapse = " "),"clusters",sep = " "))}




  #----------------------   Initialisation of the algorithm

  niter.search <- 0;
  ICL.c <- -Inf;
  classif.c <- classif.init;
  classif.new <- classif.c
  estim.new <- data$estime(classif.init);


  param.new <- estim.new$param_estim
  classif.new <- lapply(1:data$Q,function(q){Z_q = max.col(param.new$tau[[q]])
  Z_q = match(Z_q, unique(sort(Z_q)))
  names(Z_q) <- data$names_ind[[q]]; return(Z_q)})

  estim.new$param_estim$Z = classif.new
  ICL.new <- estim.new$ICL

  if (verbose) {
    mess <- paste(round(c(calc_vK(classif.new))),collapse = " " )
    print(paste("ICL :",round(ICL.new,2),". Nb of clusters: ", mess,sep = " "))
  }

  vec.ICL = ICL.new
  RES = list()
  RES[[niter.search + 1]] <- estim.new;
  #----------------------------------------------  ALGORITHM
  while((ICL.new > ICL.c)&(niter.search<1000)){
    #
    niter.search <- niter.search + 1;

    #
    ICL.c <- ICL.new
    classif.c <- classif.new

    # list of new classif deriving from the plitting of one cluster or the merging of 2 clusters in cluterisation classif.c
    # (These clusterisations will serve as initisalition of the EM algorithm for models)
    list_classif_init <- sequential_initialize(classif.c ,data,Kmin,Kmax,os);
    L = length(list_classif_init)

    #vK.c = calc_vK(classif.c)


   #  diff_vK = sapply(list_classif_init,function(cl){
   #    return(sum(abs(calc_vK(cl) - vK.c)))
   #  })
   #
   #  # remove bad initializations
   #
   #  ibadclus = which(diff_vK!=1)
   #  if (length(ibadclus)>0)
   #  {
   #    list_classif_init = list_classif_init[-ibadclus]
   #    L = length(list_classif_init)
   #  }
   #
   #  print(ibadclus)
   # browser()

    if (os == "Windows") {
      all_estim <- lapply(1:L,function(l){
        #print(l);
        estim.c.l <- data$estime(list_classif_init[[l]])})
    }else{
      all_estim <- mclapply(1:L,function(l){estim.c.l <- data$estime(list_classif_init[[l]])},mc.cores = nb_cores)
    }



    all_estim <- data$clean_results(all_estim)
    ICL.new <- all_estim[[1]]$ICL
    vec.ICL <- c(vec.ICL,ICL.new)

    if (ICL.new > ICL.c) {

      estim.new <- all_estim[[1]]
      estim.new$param_estim$Z <- lapply(1:data$Q,function(q){Z_q <- max.col(estim.new$param_estim$tau[[q]]);
      Z_q = match(Z_q, unique(sort(Z_q)))
      names(Z_q) <- data$names_ind[[q]]; return(Z_q)})

      RES[[niter.search + 1]] <- estim.new;
      param.new <- estim.new$param_estim
      classif.new <- param.new$Z
      if (verbose) {
        mess <- paste(round(c(calc_vK(classif.new))),collapse = " " )
        print(paste("ICL :",round(ICL.new,2),". Nb of clusters: ", mess,sep = " "))
      }
    }

  }

  return(RES)
}

###########################"




