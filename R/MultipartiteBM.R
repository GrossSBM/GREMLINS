#' Model selection and estimation of MBM
#'
#' Select the number of blocks per functional group using a stepwise search and estimate parameters
#'
#' @param listNet A list of networks (defined via the function DefineNetwork) i.e. multipartite networks
#' @param namesfg Names of functional groups (FG) (must correspond to names in listNet)
#' @param vKmin A vector of minimal number of blocks per functional group provided in the same order as in namesfg.
#'              vKmin can be a single value (same minimal number of blocks for all the FGs) or a vector with size equal to the number of FGs
#'              If vKmin is not specified,  vKmin = 1.
#' @param vKmax A vector of maximal number of blocks per functional group provided in the same order as in namesfg.
#'              vKmin can be a single value (same minimal number of blocks for all the FGs) or a vector with size equal to the number of FGs
#'              If vKmin is not specified,  vKmin = 1.
#' @param vKinit A vector of initial number of blocks per functional group provided in the same order as in namesfg.
#'               if vKinit is not specified, then several initialisations will be used :  vKinit = vKmin, and vKinit = floor((vKmax + vKmin)/2)
#' @param init.BM If init.BM   =  TRUE, then an aditional initialisation is done using simple LBM or SBM on each network separatly. The default value is FALSE
#' @param save Set to TRUE to save the estimated parameters for intermediate visited models. Otherwise, only the better model (in ICL sense) is the ouput
#' @param verbose Set to TRUE to display the current step of the search algorithm
#' @return a list of estimated parameters for the different models ordered by decreasing ICL. If save=FALSE, the length is of length 1
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
#' list_networks <- rMBM(v_NQ ,E , type_inter, vdistrib, lpi, ltheta, seed=NULL, namesfg= c('A','B','D'))
#' res <- MultipartiteBM(list_networks,namesfg = NULL, vdistrib = c('bernoulli",'poisson','poisson'), vKmin = 1,vKmax = 10,vKinit = NULL,verbose = TRUE, save=FALSE)
#' res2 <- MultipartiteBM(list(Agr,Bgr),namesfg = c("1","2"),vKmin = c(1,1),vKmax = c(10,10),vKinit = NULL,init.BM = TRUE, save=FALSE, verbose = TRUE)
#' @export

MultipartiteBM = function(listNet, namesfg = NULL, vdistrib = NULL , vKmin = 1 , vKmax = 10 , vKinit = NULL , init.BM = FALSE , save=FALSE , verbose = TRUE)
{


  #----------------- Formatting the data ---
  dataR6 = FormattingData(listNet,vdistrib)

  #------------------------------------------
  if (verbose) {
    Nb.entities <- dataR6$v_NQ;
    names(Nb.entities) <- dataR6$namesfg;
    print("------------Nb of entities in each functional group--------------")
    print(Nb.entities)
  }

  #------------------- Check the order of names_FG
  if (dataR6$Q == 1) {
    namesfg <- dataR6$namefg
  }else {# dataR6$Q > 1
    if ( (length(vKmin) == dataR6$Q) & is.null(namesfg))  {stop("Please specify the names of the Functional Groups")}
    if ( (length(vKmax) == dataR6$Q) & is.null(namesfg))  {stop("Please specify the names of the Functional Groups")}
    if ( (length(vKinit) == dataR6$Q) & is.null(namesfg)) {stop("Please specify the names of the Functional Groups")}
  }

  if (!is.null(namesfg)  &  !setequal(namesfg,dataR6$namesfg)) {stop("Unmatching names of Functional Groups")}

  if (is.null(vKmin)) {vKmin <- 1; print("The minimum number of clusters has been set to 1")}
  if (is.null(vKmax)) {vKmax <- 10; print("The maximum number of clusters has been set to 10")}
  if (length(vKmin) == 1) {vKmin <- rep(vKmin,dataR6$Q)}else{if (length(vKmin) != dataR6$Q) {stop("Lower bounds on vK are not of the adequate size")}}
  if (length(vKmax) == 1) {vKmax <- rep(vKmax,dataR6$Q)}else{if (length(vKmax) != dataR6$Q) {stop("Upper bounds on vK are not of the adequate size")}}

  #------------------- Reorder the vKmin, vKmax and vKinit to match the order of  dataR6$namesfg
  vKmin_permut <- vKmin
  vKmax_permut <- vKmax
  vKinit_permut <- vKinit
  for (q in 1:dataR6$Q) {
    wq <- which(dataR6$namesfg == namesfg[q])
    vKmin_permut[wq] <- vKmin[q]
    vKmax_permut[wq] <- vKmax[q]
    if ((length(vKinit) == dataR6$Q) &  length(vKinit) > 1) {vKinit_permut[q] <- vKinit[wq]}
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
    if (vKmax[q] > dataR6$v_NQ[q])
    {
      vKmax[q] = dataR6$v_NQ[q]
      print(paste("Kmax[",q,"]  was set to ",dataR6$v_NQ[q],sep = ""))
    }
  }



  if (is.null(vKinit)) {
    vKinit_list <- list(vKmin)
    vKmean <- floor((vKmax + vKmin)/2)
    if (any(vKmean != vKmin)) { vKinit_list[[2]] <- vKmean }
  }else{vKinit_list <- list(vKinit)}




  ################ ESTIMATION starting from one or two initialisation
  # classif CAH
  param.init <- genBMfit$new(vK = vKinit_list[[1]],vdistrib = dataR6$vdistrib)
  classif.init = initialize(dataR6,param.init,method = "CAH")$groups
  R = dataR6$search_nb_clusters(classif.init,Kmin = vKmin,Kmax = vKmax,verbose = verbose)
  if (length(vKinit_list) > 1) {
    param.init <- genBMfit$new(vK = vKinit_list[[2]],vdistrib = dataR6$vdistrib)
    classif.init = initialize(dataR6,param.init,method = "CAH")$groups
    R <- c(R,dataR6$search_nb_clusters(classif.init,Kmin = vKmin,Kmax = vKmax,verbose = verbose))}

  if (init.BM)
  {
    if (dataR6$card_E == 1) {print("initialisation based on each network is not relevant")}
    else {
    list_classif.initBM = lapply(1:dataR6$Q,function(q){list()})
    names(list_classif.initBM) = dataR6$namesfg

    lapply(1:dataR6$card_E, function(e){
      #version GREMLIN
      if (dataR6$type_inter[e] == "inc") { indFG = dataR6$E[e,]} else {indFG = dataR6$E[e,1]}
      estim = MultipartiteBM(list(listNet[[e]]),namesfg = dataR6$namesfg[indFG] , vKmin = vKmin[indFG] ,vKmax = vKmax[indFG] ,vKinit = vKmin[indFG], verbose = FALSE)
      if (dataR6$type_inter[e] == "inc")
      {
        list_classif.initBM[[dataR6$E[e,1]]] <<- c(list_classif.initBM[[dataR6$E[e,1]]],list(estim$fitted.model[[1]]$param_estim$Z[[1]]))
        list_classif.initBM[[dataR6$E[e,2]]] <<- c(list_classif.initBM[[dataR6$E[e,2]]],list(estim$fitted.model[[1]]$param_estim$Z[[2]]))
      } else {
        list_classif.initBM[[dataR6$E[e,1]]] <<- c(list_classif.initBM[[dataR6$E[e,1]]],list(estim$fitted.model[[1]]$param_estim$Z[[1]]))
      }
      #version blockmodels
    #   estim =  switch(dataR6$type_inter[e],
    #          "inc" = BM_bernoulli("LBM",dataR6$mats[[e]],verbosity=0,plotting=""),
    #          "diradj" = BM_bernoulli("SBM",dataR6$mats[[e]],verbosity=0,plotting=""),
    #          "adj" = BM_bernoulli("SBM_sym",dataR6$mats[[e]],verbosity=0,plotting=""))
    #
    #   estim$estimate()
    #   k = which.max(estim$ICL)
    #   best_clust = estim$memberships[[k]]
    #
    #   if (dataR6$type_inter[e]=="inc")
    #   {
    #       list_classif.initBM[[dataR6$E[e,1]]] <<- c(list_classif.initBM[[dataR6$E[e,1]]],list(apply(best_clust$Z1,1,which.max)))
    #       list_classif.initBM[[dataR6$E[e,2]]] <<- c(list_classif.initBM[[dataR6$E[e,2]]],list(apply(best_clust$Z2,1,which.max)))
    #   } else {
    #     list_classif.initBM[[dataR6$E[e,1]]] <<- c(list_classif.initBM[[dataR6$E[e,1]]],list(apply(best_clust$Z,1,which.max)))
    #   }
     })



    Nb_classif.initBM = lapply(list_classif.initBM,function(l) 1:length(l))
    combin_classif.initBM = as.matrix(expand.grid(Nb_classif.initBM))


    lapply(1:nrow(combin_classif.initBM),function(i)
    {
         rowcombin = as.vector(combin_classif.initBM[i,])
         classif.init = lapply(1:dataR6$Q, function(q) list_classif.initBM[[q]][[rowcombin[q]]])
          R <<- c(R,dataR6$search_nb_clusters(classif.init,Kmin = vKmin,Kmax = vKmax,verbose = verbose))
    }
   )
  }
  }




  #-------------------- cleaning the results
  res <- dataR6$clean_results(R) # remove models that have been estimated twice or more to keep the estimation with the better J

  lapply(1:length(res),function(k){names(res[[k]]$param_estim$vK) <<- namesfg})



  ############### PRINT RESULT#######################
  if (verbose) {
    mess <- paste(res[[1]]$param_estim$vK,collapse = " " )
    mess <- paste("Best model------ ICL :",round(res[[1]]$ICL,2),". Nb of clusters: ", mess,sep = " ")
    mess <- paste(mess, "for",paste(dataR6$namesfg,collapse = " , " ),"respectively",sep = ' ')
    print(mess)
  }

  ############# RESULTATS #############################"
  if (save) {return(list(fitted.model = res ,listNet = listNet))}else{return(list(fitted.model = list(res[[1]]) ,listNet = listNet))}

}

