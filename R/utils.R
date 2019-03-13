# ------------------ Test if a number is an integer
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {abs(x - round(x)) < tol}


# ------------------ Test if a number is a positive integer
is.poswholenumber <- function(x, tol = .Machine$double.eps^0.5)  {is.wholenumber(x) & (x >= 0)}


#---------------------- ARI
adjustedRandIndex <- function (x, y)
{
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y))
    stop("arguments must be vectors of the same length")
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1)))
    return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b +
      a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}

# ------------------ Laplace distribution
rlaplace <- function(n, location = 0, scale = 1){
  mu <- c(location)
  b <- c(scale)
  U <- runif(n,-0.5,0.5)
  res <- mu - b * sign(U) * log(1 - 2 * abs(U))
  return(res)
}

dlaplace <- function(x, location = 0, scale = 1, log = FALSE){
  mu <- c(location)
  b <- c(scale)
  log.df <- -log(2 * b) - abs(x - mu) / b;
  if (log == TRUE) {return(log.df)} else {return(exp(log.df))}

}

argminWeightedTAV <- function(Y,weights,o = NULL,isYordered = FALSE){

    if (length(Y) != length(weights) ){stop('pb of vector size in argminWeightedTAV')}
    W <- weights
    if (is.null(o )) {o <- order(Y,decreasing = FALSE)}
    if (!isYordered) {Yo <- Y[o]} else {Yo <- Y}
    A <- sum(cumsum(W[o]) < sum(W[o]) * 0.5)
    return(Yo[A + 1])
  }


################################################################################################################
###----------------Check the sizes of the matrix with respect to the Functional group they involve  -------#####
################################################################################################################

checkExtract = function(list_Mat,matE)
  #same inputs as in VEM_gen_BM
{

  #number of functional groups
  functionalGroups = unique(as.vector(matE))
  functionalGroups = functionalGroups[functionalGroups > 0]

  #extracting dim of matrices


  where = lapply(functionalGroups,function(q){
    a = which(matE == q,arr.ind = TRUE)
  })

  #check SBM sym
  indSBM <- which(matE[,2] == 0)
  vsym <- sapply(indSBM,function(i){!isSymmetric(list_Mat[[i]])})
  if (length(vsym) > 0) {
    if (sum(vsym) > 0)
    {
      stop("symmetry problem with one ore more of the interaction matrix")
    }
  }

  #check SBM nonsym
  indSBM = which(matE[,2] == -1)
  vnonsym = sapply(indSBM,function(i){
    calc = dim(list_Mat[[i]])
    return(calc[1] != calc[[2]])
  })
  if (length(vnonsym) > 0) {
    if (sum(vnonsym) > 0)
    {
      stop("one ore more of intra interaction matrices is not square")
    }
  }


  #extract n_q
  n_qs = lapply(where,function(w_q)
  {
    sapply(list_Mat[w_q[,1]],dim)[cbind(w_q[,2],1:nrow(w_q))]
  })


  #check n_q
  n_q = sapply(n_qs,function(n)
  {
    temp = unique(n)
    if (length(temp) != 1)
    {stop("dimensions mismatch")}
    return((unique(n)))
  })
  #return
  return(list(n_q = n_q,where_q = where))
}

transfoE <- function(E, typeInter){
  Ecode <- E
  Ecode[which(typeInter == "adj"),2] <- 0
  Ecode[which(typeInter == "diradj"),2] <- -1
  return(Ecode)
}

################################################################################################################
#------------------ transform the data given by the user into  a R6 type coll_interaction object
################################################################################################################

formattingData <- function(list_Net,v_distrib = NULL)
{


  mats = lapply(list_Net,function(net){return(net$mat)})
  intFG = lapply(list_Net,function(net){return(c(net$rowFG,net$colFG))})
  intFG = do.call(rbind,intFG)
  types = sapply(list_Net,function(net){return(net$typeInter)})
  return(CollInteraction$new(mats = mats,E_FG = intFG,typeInter = types,v_distrib = v_distrib))
}


##############################################################################################################
###-----------------------------------FUNCTIONS USED in the VEM algortihm ---------------
##############################################################################################################

#----------------------  distance on tau

distTau  <- function(tau,tauOld)
{
  Q <- length(tau)
  vdis <- sapply(1:Q,function(q){
    return(sqrt(sum(as.vector(tau[[q]] - tauOld[[q]])^2)))
  })
  return(sum(vdis))
}


#------------------------------ prevent the parameters from being too close from the bound of the domain


readjustPi <- function(pi,eps)
{
  pi[pi < eps] <- eps
  pi[pi > (1 - eps)] <- 1 - eps
  return(pi)
}

readjustTheta <- function(theta,eps, distrib)
{
  if (distrib == 'bernoulli') {
    theta[theta < eps] = eps
    theta[theta > 1 - eps] = 1 - eps }
  if (distrib == 'poisson') { theta[theta < eps] = eps }
  if (distrib == 'gaussian'){ theta$sd[theta$sd < eps] = eps }
  if (distrib == 'laplace') { theta$scale[theta$scale < eps] = eps }
  return(theta)
}


#------------------- computing ICL and likelihood
compLikICLInt = function(tau,list_theta,list_pi,matE,list_Mat,n_q,v_K,v_distrib)
{

  #browser()
  cardE <- nrow(matE)
  Q <-  length(list_pi)

  #condLik
  condLiks = sapply(1:cardE,function(e)
  {
    gr = matE[e,1]
    gc = matE[e,2]
    don = list_Mat[[e]]
    if (v_distrib[e] == 'bernoulli') {Unmdon = 1 - don}
    facteur = 1
    if (gc < 1)
    {
      if (gc == 0)  facteur = 1/2 #sbm sym
      gc = gr
      if (v_distrib[e] == 'bernoulli') {diag(Unmdon) = 0}
    }
    if (v_distrib[e] == 'bernoulli') {
      prov = (tau[[gr]]) %*% log(list_theta[[e]]) %*% t(tau[[gc]])
      prov1m = (tau[[gr]]) %*% log(1 - list_theta[[e]]) %*% t(tau[[gc]])
      return((sum(don * prov) + sum(Unmdon * prov1m)) * facteur)
    }
    if (v_distrib[e] == 'poisson') {
      prov = (tau[[gr]]) %*% log(list_theta[[e]]) %*% t(tau[[gc]])
      prov2 = (tau[[gr]]) %*%  list_theta[[e]]  %*% t(tau[[gc]])
      Unit <- matrix(1,nrow = nrow(don),ncol = ncol(don))
      return((sum(don * prov) - sum(Unit * prov2)) * facteur)
    }

  }
  )
  condLik = sum(condLiks)

  margLiks = sapply(1:Q,function(q)
  {
    return(sum(tau[[q]]*matrix(log(list_pi[[q]]),nrow(tau[[q]]),ncol(tau[[q]]),byrow = TRUE)))
  })
  margLik = sum(margLiks)

  entros = sapply(1:Q,function(q){return(-sum(log(tau[[q]]) * tau[[q]]))})
  entro = sum(entros)

  #penalty
  penMats = vapply(1:cardE,function(s)
  {
    gr = matE[s,1]
    gc = matE[s,2]
    if (gc > 0) #LBM
    {
      return(c(v_K[gr] * v_K[gc], prod(dim(list_Mat[[s]]))))
    }
    if (gc == 0) #SBM sym
    {
      return(c(v_K[gr] * (v_K[gr] + 1) / 2, nrow(list_Mat[[s]]) * (nrow(list_Mat[[s]]) - 1)/2))
    }
    if (gc == -1) #SBM non sym
    {
      return(c(v_K[gr] * (v_K[gr]), nrow(list_Mat[[s]]) * (nrow(list_Mat[[s]]) - 1)))
    }

  },c(1.0,2.0)
  )


  penICL = sum((v_K - 1) * log(n_q))  + sum(penMats[1,]) * log(sum(penMats[2,]))


  return(list(condLik = condLik,margLik = margLik,entr = entro,pen = penICL))

}


#distance between pis
# distListPi = function(list_pi,list_piOld)
# {
#   Q <- length(list_pi)
#   v_dis <-  sapply(1:Q,function(q){
#     return(sqrt(sum(as.vector(list_pi[[q]] - list_piOld[[q]])^2)))
#   })
#   return(sum(v_dis))
# }

#distance between l_theta
distListTheta <- function(list_theta,list_thetaOld)
{
  Q <- length(list_theta)
  v_dis <- sapply(1:Q,function(q){
    return(sqrt(sum(as.vector(list_theta[[q]] - list_thetaOld[[q]])^2)))
  })
  return(sum(v_dis))
}






#verif
# ARIS <- function(res,classifVrai)
# {
#   Q <- length(res$tau)
#   a <- sapply(1:Q,function(q){
#     cl=apply(res$tau[[q]],1,which.max)
#     adjustedRandIndex(cl,classifVrai[[q]])})
#   return(a)
# }
#
# # frobenius <- function(X){
#   return(sqrt(sum((X)^2)))
# }

# erparam <- function(res,paramVrai)
# {
#   Q <- length(res$alpha)
#   cardE <- length(res$pi)
#
#   eralpha=sapply(1:Q,function(q){
#     C=paramVrai$alpha[[q]]
#     alph=res$alpha[[q]]
#     K_q=length(res$alpha[[q]])
#     o=min(sapply(permn(K_q),
#                  function(perm) {
#                    A <- alph[perm]
#                    return(frobenius(C-A)/frobenius(C))}))
#     return(o)
#   })
#
#   erpi=sapply(1:cardE,function(q)
#   {
#     C=paramVrai$pi[[q]]
#     pi=res$pi[[q]]
#     K_q=nrow(res$pi[[q]])
#     K_qprime=ncol(res$pi[[q]])
#     o=min(sapply(permn(K_q),
#                  function(perm) {
#                    oprime=min(sapply(permn(K_qprime),function(perm2){
#                      A <- pi[perm,perm2]
#                      return(frobenius(C-A)/frobenius(C))}))
#
#                    return(oprime)})
#     )
#
#     return(o)
#
#   })
#
#   return(list(eralpha,erpi))
# }






##############################################################################################################
# ------------------ Cleaning a list of estimations
##############################################################################################################

cleanEstim <- function(dataR6,R){


  ### first  : for equal model keep the estimate wit the highest J
  seq_J <- sapply(R ,function(u){max(u$vJ)})
  o <- order(seq_J,decreasing = TRUE)
  R.ordered.1 <- lapply(o,function(i){R[[i]]})

  if (dataR6$Q == 1) {
  seq_NbClust <- cbind((sapply(R.ordered.1,function(u){u$paramEstim$v_K})),seq_J[o],1:length(R))
  } else {
  seq_NbClust <- cbind(t(sapply(R.ordered.1,function(u){u$paramEstim$v_K})),seq_J[o],1:length(R))
  }


  if (length(R) > 1)
  {
    seq_NbClust <- seq_NbClust[!duplicated(seq_NbClust[,1:dataR6$Q]),]
    if (is.null(nrow(seq_NbClust))) {seq_NbClust = matrix(seq_NbClust,nrow = 1)}
    R <-  R.ordered.1[seq_NbClust[,dataR6$Q + 2]]
  }

  ### Order the results by ICL
  seq_ICL <- sapply(R,function(u){u$ICL})
  o <- order(seq_ICL,decreasing = TRUE)
  res <- lapply(o,function(i){R[[i]]})
  return(res)
}


##############################################################################################################
# ------------------ Cleaning a list of initialisations
##############################################################################################################


#---------- COMPARISON of two classfications.

testClassif <- function(classif1,classif2){

  if (length(classif1) != length(classif2) )
  {
    stop('classification can not be compared because objects of different sizes')
  }
  ARI <- mean(vapply(1:length(classif1),function(q){adjustedRandIndex(classif1[[q]],classif2[[q]])},1))
  if (ARI == 1) {return(TRUE)} else {return(FALSE)}
}


# ----------------- Compare the initialisation on Z where we already started from

cleanCollectionClassif <- function(collectionClassif,indRef){

  L <- length(collectionClassif)
  matTest <- matrix(NA,L,L)
  for (l in 2:L) {for (k in 1:(l - 1)) {matTest[l,k] <- testClassif(collectionClassif[[k]],collectionClassif[[l]])}}
  w <- which(rowSums(matTest[indRef:L,],na.rm = TRUE) > 0) - indRef +  1
  u <- (1:L)[-w]
  res <- lapply(u,function(i){collectionClassif[[i]]})
  return(res)
}










