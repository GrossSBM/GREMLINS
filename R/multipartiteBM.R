#' Model selection and parameter estimation of MBM
#'
#' Select the number of blocks and identify the blocks per functional group using a variational EM algorithm
#'
#' @param list_Net a list of networks (defined via the function defineNetwork) i.e. a multipartite network
#' @param v_distrib an optional vector of characters of length the number of networks and specifying the distribution used in each network (possible values \code{bernoulli,poisson,gaussian,laplace}).   If not provided, the model will be 'bernoulli' for all the interactions  matrices.
#' @param namesFG  an optional vector of  characters containing the names of functional groups (FG) (If Specified, must correspond to the  names in \code{list_Net}).
#' @param v_Kmin an optional vector of integers, specifying the minimal number of blocks per functional group (must be provided in the same order as in \code{namesFG}).
#'              v_Kmin may be a single value (same minimal number of blocks for all the FGs) or a vector with size equal to the number of FGs.
#'              Default value  \code{= 1}.
#' @param v_Kmax an optional vector of integers specifying the maximal number of blocks per functional group provided in the same order as in \code{namesFG}.
#'              v_Kmax may be a single value (same maximal number of blocks for all the FGs) or a vector with size equal to the number of FGs.
#'              Default value  \code{= 10}.
#' @param v_Kinit an optional vector of integers specifying initial numbers of blocks per FG provided in the same order as in \code{namesFG}.
#'               if \code{v_Kinit} is not specified, then   \code{v_Kinit = v_Kmin}
#' @param initBM an optional boolean. If initBM = TRUE  an aditional initialisation is done using simple LBM or SBM on each network separatly.
#'               Default value  \code{= TRUE}
#' @param keep an optional boolean. If TRUE  return the estimated parameters for intermediate visited models. Otherwise, only the better model (in ICL sense) is the ouput. Default value \code{= FALSE}.
#' @param verbose an optional boolean. If  TRUE, display the current step of the search algorithm
#' @param nbCores an optional integer specifying the number or cores used for the estimation. Not parallelized on windows. If \code{ncores = NULL}, then half of the cores are used.
#' @param maxiterVE an optional integer  specifying the maximum number of iterations in the VE step of the VEM algorithm. If NULL then default value  \code{= 1000}
#' @param maxiterVEM an optional integer  specifying the maximum number of iterations of the VEM algorithm. If NULL then default value Default value  \code{= 1000}
#' @details The function \code{multipartiteBM} selects the better numbers of blocks in each FG (with a penalized likelihood criterion). The model selection is performed with a forward backward strategy and the likelihood of each model is maximized with a variational EM).

#'
#' @return a list of estimated parameters for the different models ordered by decreasing ICL. If keep \code{= FALSE}, the length is of length 1 (only the better model is returned).
#' \describe{
#'   \item{\code{fittedModel}}{contains the results of the inference. \code{res$fittedModel[[1]]}  is a list with fields
#'   \describe{
#'   \item{\code{paramEstim}}{a MBMfit object.}
#'   \item{\code{ICL}}{the penalized likelihood criterion ICL.}
#'   \item{\code{vJ}}{the sequence of the varational bound of the likelihood through iterations of the VEM.}
#'   \item{\code{convergence}}{TRUE if the VEM reached convergence.}
#'   }
#'   }
#'   \item{\code{list_Net}}{ contains the data.}
#'}
#'
#'
#'
#'
#'
#'
#' @examples

#' namesFG <- c('A','B')
#' list_pi <- list(c(0.5,0.5),c(0.3,0.7)) # prop of blocks in each FG
#' E  <-  rbind(c(1,2),c(2,2)) # architecture of the multipartite net.
#' typeInter <- c( "inc","diradj")
#' v_distrib <- c('gaussian','bernoulli')
#' list_theta <- list()
#' list_theta[[1]] <- list()
#' list_theta[[1]]$mean  <- matrix(c(6.1, 8.9, 6.6, 3), 2, 2)
#' list_theta[[1]]$var  <-  matrix(c(1.6, 1.6, 1.8, 1.5),2, 2)
#' list_theta[[2]] <- matrix(c(0.7,1.0, 0.4, 0.6),2, 2)
#' list_Net <- rMBM(v_NQ = c(30,30),E , typeInter, v_distrib, list_pi,
#'                 list_theta, namesFG = namesFG, seed = 2)$list_Net
#' res_MBMsimu <- multipartiteBM(list_Net, v_distrib,
#'                               namesFG = c('A','B'), v_Kinit = c(2,2),
#'                               nbCores = 2,initBM = FALSE)

#' @export

multipartiteBM = function(list_Net,  v_distrib = NULL ,namesFG = NULL, v_Kmin = 1 , v_Kmax = 10 , v_Kinit = NULL , initBM = TRUE , keep = FALSE , verbose = TRUE, nbCores = NULL, maxiterVE = NULL , maxiterVEM = NULL)
{


  if ( all(v_Kmin == v_Kmax)) {stop('v_Kmin = v_Kmax. Use the function "multipartiteBMFixedModel" instead')}


  #----------------- Formatting the data ---
  dataR6 = formattingData(list_Net,v_distrib)

  if (verbose) {
    NBEntities <- dataR6$v_NQ;
    names(NBEntities) <- dataR6$namesFG;
    print("------------Nb of entities in each functional group--------------")
    print(NBEntities)


    print("------------Probability distributions on each network--------------")
    print(v_distrib)
    print("-------------------------------------------------------------------")
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

  if (dataR6$Q >1) {
    if (!is.null(v_Kinit) &  (length(v_Kinit) != dataR6$Q)) {
      print("v_Kinit was not of the adequate size. The given value has not been taken into account");
      print("-------------------------------------------------------------------")
      v_Kinit = NULL}
  }

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
    #v_Kmean <- floor((v_Kmax + v_Kmin)/2)
    #if (any(v_Kmean != v_Kmin)) { v_Kinit_list[[2]] <- v_Kmean }
  }else{v_Kinit_list <- list(v_Kinit)}




  #----------------------   ESTIMATION starting from one (given) or two initialisations  (v_Kmean and v_Kmin)

  collectionTestedClassifInit <- list()

  paramInit <- MBMfit$new(v_K = v_Kinit_list[[1]],v_distrib = dataR6$v_distrib)
  classifInit <- initialize(dataR6,paramInit,method = "CAH")$groups

  R <- dataR6$searchNbClusters(classifInit,Kmin = v_Kmin,Kmax = v_Kmax,pastICL = c(),verbose = verbose,nbCores = nbCores, maxiterVE = maxiterVE , maxiterVEM = maxiterVEM)
  indInit <- 1
  collectionTestedClassifInit[[indInit]] <-  classifInit
  pastICL <- sapply(R,function(e){e$ICL})



  # vkinit_list[[2]] and further  :  init classif CAH + searching from that point
  if (length(v_Kinit_list) > 1) {
    paramInit <- MBMfit$new(v_K = v_Kinit_list[[2]],v_distrib = dataR6$v_distrib)
    classifInit = initialize(dataR6,paramInit,method = "CAH")$groups
    indInit <- indInit +  1
    collectionTestedClassifInit[[indInit]] = classifInit
    R <- c(R,dataR6$searchNbClusters(classifInit,Kmin = v_Kmin,Kmax = v_Kmax,pastICL = pastICL,verbose = verbose,nbCores = nbCores, maxiterVE = maxiterVE ,  maxiterVEM = maxiterVEM))
    pastICL <- sapply(R,function(e){e$ICL})
  }

  indInit <- length(collectionTestedClassifInit)


  # Additional initialisation starting from a block model on each network
  if (initBM)
  {
    if (dataR6$cardE == 1) {print("initialisation based on each network is not relevant")}
    else {
      list_classifInitBM = lapply(1:dataR6$Q,function(q){list()})
      names(list_classifInitBM) = dataR6$namesFG

      for (e in 1:dataR6$cardE){
        if (dataR6$typeInter[e] == "inc") { indFG = dataR6$E[e,]} else {indFG = dataR6$E[e,1]}
        #------------ esim SBM ou LSB on one network
        estim = multipartiteBM(list(list_Net[[e]]),namesFG = dataR6$namesFG[indFG] ,  v_distrib = v_distrib[e], v_Kmin = v_Kmin[indFG] ,v_Kmax = v_Kmax[indFG] ,v_Kinit = v_Kmin[indFG],  initBM = FALSE, verbose = FALSE,  nbCores = nbCores, maxiterVE = maxiterVE , maxiterVEM = maxiterVEM)
        if (dataR6$typeInter[e] == "inc")
        {
          list_classifInitBM[[dataR6$E[e,1]]] <- c(list_classifInitBM[[dataR6$E[e,1]]],list(estim$fittedModel[[1]]$paramEstim$Z[[1]]))
          list_classifInitBM[[dataR6$E[e,2]]] <- c(list_classifInitBM[[dataR6$E[e,2]]],list(estim$fittedModel[[1]]$paramEstim$Z[[2]]))
        } else {
          list_classifInitBM[[dataR6$E[e,1]]] <- c(list_classifInitBM[[dataR6$E[e,1]]],list(estim$fittedModel[[1]]$paramEstim$Z[[1]]))
      }
     }





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
        #browser()
        R1 <- dataR6$searchNbClusters(collectionTestedClassifInit[[i]],Kmin = v_Kmin,Kmax = v_Kmax,pastICL = pastICL,verbose = verbose,nbCores = nbCores, maxiterVE = maxiterVE ,   maxiterVEM = maxiterVEM)
        R <- c(R,R1)
        pastICL <- c(pastICL,sapply(R1,function(e){e$ICL}))
       }
    }

  }
  }




  #-------------------- cleaning the results


  res <- dataR6$cleanResults(R) # remove models that have been estimated twice or more to keep the estimation with the better J

  for (k in 1:length(res)){names(res[[k]]$paramEstim$v_K) <- namesFG}




  ############### PRINT RESULT#######################
  if (verbose) {
    mess <- paste(res[[1]]$paramEstim$v_K,collapse = " " )
    mess <- paste("Best model------ ICL :",round(res[[1]]$ICL,2),". Nb of clusters: [", mess,"]",sep = " ")
    mess <- paste(mess, "for [",paste(dataR6$namesFG,collapse = " , " ),"] respectively",sep = ' ')
    print(mess)
  }

  ############# RESULTATS #############################"
  if (keep) {return(list(fittedModel = res ,list_Net = list_Net))}else{return(list(fittedModel = list(res[[1]]) ,list_Net = list_Net))}

}

