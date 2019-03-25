#' Model selection and parameter estimation of MBM
#'
#' Select the number of blocks per functional group using a stepwise search and estimate parameters
#'
#' @param list_Net A list of networks (defined via the function defineNetwork) i.e. multipartite networks
#' @param namesFG Names of functional groups (FG) (must correspond to names in list_Net)
#' @param v_distrib  Type of proababilistic distributions in each network : if 0/1 then Bernoulli, if counting then Poisson. My default  = Bernoulli.
#'                   Must give a vector wos length the number of networks in list_Net
#' @param v_Kmin A vector of minimal number of blocks per functional group provided in the same order as in namesFG.
#'              v_Kmin can be a single value (same minimal number of blocks for all the FGs) or a vector with size equal to the number of FGs
#'              If v_Kmin is not specified,  v_Kmin = 1.
#' @param v_Kmax A vector of maximal number of blocks per functional group provided in the same order as in namesFG.
#'              v_Kmin can be a single value (same minimal number of blocks for all the FGs) or a vector with size equal to the number of FGs
#'              If v_Kmin is not specified,  v_Kmin = 1.
#' @param v_Kinit A vector of initial number of blocks per functional group provided in the same order as in namesFG.
#'               if v_Kinit is not specified, then several initialisations will be used :  v_Kinit = v_Kmin, and v_Kinit = floor((v_Kmax + v_Kmin)/2)
#' @param initBM If initBM   =  TRUE, then an aditional initialisation is done using simple LBM or SBM on each network separatly. The default value is FALSE
#' @param save Set to TRUE to save the estimated parameters for intermediate visited models. Otherwise, only the better model (in ICL sense) is the ouput
#' @param verbose Set to TRUE to display the current step of the search algorithm
#' @param nbCores Number or cores used for the estimation. Not parallelized on windows. By default : half of the cores
#' @param maxiterVE  Maximum number of iterations if the VE step of the VEM algorithm. By default  = 100
#' @return a list of estimated parameters for the different models ordered by decreasing ICL. If save=FALSE, the length is of length 1
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
#' dataSim <-  rMBM(v_NQ ,E , typeInter, v_distrib, list_pi, list_theta, seed=NULL, namesFG= c('A','B','D'),keepClassif = FALSE)
#' list_Net <- dataSim$list_Net
#' res <- multipartiteBM(list_Net,namesFG = NULL, v_distrib = c("bernoulli","poisson","poisson"), v_Kmin = 1,v_Kmax = 10,v_Kinit = NULL,verbose = TRUE, save=FALSE, maxiterVE = NULL)
#' @export

multipartiteBM = function(list_Net, namesFG = NULL, v_distrib = NULL , v_Kmin = 1 , v_Kmax = 10 , v_Kinit = NULL , initBM = FALSE , save=FALSE , verbose = TRUE,nbCores = NULL, maxiterVE = NULL)
{


  #----------------- Formatting the data ---
  dataR6 = formattingData(list_Net,v_distrib)

  #------------------- Messages  ----------
  if (verbose) {
    NBEntities <- dataR6$v_NQ;
    names(NBEntities) <- dataR6$namesFG;
    print("------------Nb of entities in each functional group--------------")
    print(NBEntities)


    print("------------Probability distributions on each network--------------")
    print(v_distrib)
  }

  #------------------- Check the order of names_FG
  if (dataR6$Q == 1) {
    namesFG <- dataR6$nameFG
  }else {# dataR6$Q > 1
    if ( (length(v_Kmin) == dataR6$Q) & is.null(namesFG))  {stop("Please specify the names of the Functional Groups")}
    if ( (length(v_Kmax) == dataR6$Q) & is.null(namesFG))  {stop("Please specify the names of the Functional Groups")}
    if ( (length(v_Kinit) == dataR6$Q) & is.null(namesFG)) {stop("Please specify the names of the Functional Groups")}
  }

  if (!is.null(namesFG)  &  !setequal(namesFG,dataR6$namesFG)) {stop("Unmatching names of Functional Groups")}

  if (is.null(v_Kmin)) {v_Kmin <- 1; print("The minimum number of clusters has been set to 1")}
  if (is.null(v_Kmax)) {v_Kmax <- 10; print("The maximum number of clusters has been set to 10")}
  if (length(v_Kmin) == 1) {v_Kmin <- rep(v_Kmin,dataR6$Q)}else{if (length(v_Kmin) != dataR6$Q) {stop("Lower bounds on v_K are not of the adequate size")}}
  if (length(v_Kmax) == 1) {v_Kmax <- rep(v_Kmax,dataR6$Q)}else{if (length(v_Kmax) != dataR6$Q) {stop("Upper bounds on v_K are not of the adequate size")}}

  #------------------- Reorder the v_Kmin, v_Kmax and v_Kinit to match the order of  dataR6$namesFG
  v_Kmin_permut <- v_Kmin
  v_Kmax_permut <- v_Kmax
  v_Kinit_permut <- v_Kinit
  for (q in 1:dataR6$Q) {
    wq <- which(dataR6$namesFG == namesFG[q])
    v_Kmin_permut[wq] <- v_Kmin[q]
    v_Kmax_permut[wq] <- v_Kmax[q]
    if ((length(v_Kinit) == dataR6$Q) &  length(v_Kinit) > 1) {v_Kinit_permut[q] <- v_Kinit[wq]}
  }
  v_Kinit <- v_Kinit_permut
  v_Kmin <- v_Kmin_permut
  v_Kmax <- v_Kmax_permut




  #------------------------  Bounds of the numer of clusters
  if (is.null(v_Kmin)) {v_Kmin = 1; print("The minimum number of clusters has been set to 1")}
  if (is.null(v_Kmax)) {v_Kmax = 10; print("The maximum number of clusters has been set to 10")}

  if (length(v_Kmin) == 1) {v_Kmin = rep(v_Kmin,dataR6$Q)} else {if (length(v_Kmin) != dataR6$Q) {stop("Lower bounds on v_K are not of the adequate size")}}
  if (length(v_Kmax) == 1) {v_Kmax = rep(v_Kmax,dataR6$Q)} else {if (length(v_Kmax) != dataR6$Q) {stop("Upper bounds on v_K are not of the adequate size")}}

  for (q in 1:dataR6$Q)
  {
    if (v_Kmax[q] > dataR6$v_NQ[q])
    {
      v_Kmax[q] = dataR6$v_NQ[q]
      print(paste("Kmax[",q,"]  was set to ",dataR6$v_NQ[q],sep = ""))
    }
  }


  #------------------------  Initial values for the number of the  numbers  of clusters VK
  if (is.null(v_Kinit)) {
    v_Kinit_list <- list(v_Kmin)
    v_Kmean <- floor((v_Kmax + v_Kmin)/2)
    if (any(v_Kmean != v_Kmin)) { v_Kinit_list[[2]] <- v_Kmean }
  }else{v_Kinit_list <- list(v_Kinit)}




  #----------------------   ESTIMATION starting from one (given) or two initialisations  (v_Kmean and v_Kmin)

  collectionTestedClassifInit <- list()

  paramInit <- MBMfit$new(v_K = v_Kinit_list[[1]],v_distrib = dataR6$v_distrib)
  classifInit = initialize(dataR6,paramInit,method = "CAH")$groups
  R = dataR6$searchNbClusters(classifInit,Kmin = v_Kmin,Kmax = v_Kmax,verbose = verbose,nbCores = nbCores, maxiterVE = maxiterVE)
  indInit = 1
  collectionTestedClassifInit[[indInit]] <-  classifInit


  # vkinit_list[[2]] and further  :  init classif CAH + searching from that point
  if (length(v_Kinit_list) > 1) {
    paramInit <- MBMfit$new(v_K = v_Kinit_list[[2]],v_distrib = dataR6$v_distrib)
    classifInit = initialize(dataR6,paramInit,method = "CAH")$groups
    indInit <- indInit +  1
    collectionTestedClassifInit[[indInit]] = classifInit
    R <- c(R,dataR6$searchNbClusters(classifInit,Kmin = v_Kmin,Kmax = v_Kmax,verbose = verbose,nbCores = nbCores, maxiterVE = maxiterVE))
  }

  indInit <- length(collectionTestedClassifInit)


  # Additional initialisation starting from a block model on each network
  if (initBM)
  {
    if (dataR6$cardE == 1) {print("initialisation based on each network is not relevant")}
    else {
      list_classifInitBM = lapply(1:dataR6$Q,function(q){list()})
      names(list_classifInitBM) = dataR6$namesFG

      lapply(1:dataR6$cardE, function(e){
        if (dataR6$typeInter[e] == "inc") { indFG = dataR6$E[e,]} else {indFG = dataR6$E[e,1]}
        estim = multipartiteBM(list(list_Net[[e]]),namesFG = dataR6$namesFG[indFG] ,  v_distrib = v_distrib[e], v_Kmin = v_Kmin[indFG] ,v_Kmax = v_Kmax[indFG] ,v_Kinit = v_Kmin[indFG],  initBM = FALSE, verbose = FALSE, maxiterVE = maxiterVE)
        if (dataR6$typeInter[e] == "inc")
        {
          list_classifInitBM[[dataR6$E[e,1]]] <<- c(list_classifInitBM[[dataR6$E[e,1]]],list(estim$fittedModel[[1]]$paramEstim$Z[[1]]))
          list_classifInitBM[[dataR6$E[e,2]]] <<- c(list_classifInitBM[[dataR6$E[e,2]]],list(estim$fittedModel[[1]]$paramEstim$Z[[2]]))
        } else {
          list_classifInitBM[[dataR6$E[e,1]]] <<- c(list_classifInitBM[[dataR6$E[e,1]]],list(estim$fittedModel[[1]]$paramEstim$Z[[1]]))
      }
     })




      Nb_classifInitBM = lapply(list_classifInitBM,function(l) 1:length(l))
      combin_classifInitBM = as.matrix(expand.grid(Nb_classifInitBM))
      indRef <- indInit
      for (i in 1:nrow(combin_classifInitBM))
      {
        indInit <- indInit + 1
        rowcombin = as.vector(combin_classifInitBM[i,])
        classifInit = lapply(1:dataR6$Q, function(q) list_classifInitBM[[q]][[rowcombin[q]]])
        collectionTestedClassifInit[[indInit]] <- classifInit
    }

    collectionTestedClassifInit <- cleanCollectionClassif(collectionTestedClassifInit,indRef)
    if (length(collectionTestedClassifInit) > indRef) {
      for (i in (indRef + 1):length(collectionTestedClassifInit))
      {
         R <<- c(R,dataR6$searchNbClusters(classifInit,Kmin = v_Kmin,Kmax = v_Kmax,verbose = verbose,nbCores = nbCores, maxiterVE = maxiterVE))
      }
    }

  }
  }




  #-------------------- cleaning the results

  res <- dataR6$cleanResults(R) # remove models that have been estimated twice or more to keep the estimation with the better J

  lapply(1:length(res),function(k){names(res[[k]]$paramEstim$v_K) <<- namesFG})



  ############### PRINT RESULT#######################
  if (verbose) {
    mess <- paste(res[[1]]$paramEstim$v_K,collapse = " " )
    mess <- paste("Best model------ ICL :",round(res[[1]]$ICL,2),". Nb of clusters: ", mess,sep = " ")
    mess <- paste(mess, "for",paste(dataR6$namesFG,collapse = " , " ),"respectively",sep = ' ')
    print(mess)
  }

  ############# RESULTATS #############################"
  if (save) {return(list(fittedModel = res ,list_Net = list_Net))}else{return(list(fittedModel = list(res[[1]]) ,list_Net = list_Net))}

}

