# ------------------ Test if a number is an integer
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {abs(x - round(x)) < tol}


# ------------------ Test if a number is a positive integer
is.poswholenumber <- function(x, tol = .Machine$double.eps^0.5)  {is.wholenumber(x) & (x >= 0)}



################################################################################################################
###----------------Check the sizes of the matrix with respect to the Functional group they involve  -------#####
################################################################################################################

check_extract=function(list_Mat,mat_E)
  #same inputs as in VEM_gen_BM
{

  #number of functional groups
  functional_groups = unique(as.vector(mat_E))
  functional_groups = functional_groups[functional_groups>0]
  Q = length(functional_groups)

  #extracting dim of matrices
  dims = vapply(list_Mat,dim,1:2)

  where = lapply(functional_groups,function(q){
    a = which(mat_E == q,arr.ind = TRUE)
  })

  #check SBM sym
  ind_SBM=which(mat_E[,2]==0)
  vsym=sapply(ind_SBM,function(i){!isSymmetric(list_Mat[[i]])})
  if (length(vsym)>0){
    if (sum(vsym)>0)
    {
      stop("symmetry problem with one ore more of the interaction matrix")
    }
  }

  #check SBM nonsym
  ind_SBM=which(mat_E[,2]==-1)
  vnonsym=sapply(ind_SBM,function(i){
    calc=dim(list_Mat[[i]])
    return(calc[1]!=calc[[2]])
  })
  if (length(vnonsym)>0){
    if (sum(vnonsym)>0)
    {
      stop("one ore more of intra interaction matrices is not square")
    }
  }


  #extract n_q
  n_qs=lapply(where,function(w_q)
  {
    sapply(list_Mat[w_q[,1]],dim)[cbind(w_q[,2],1:nrow(w_q))]
  })


  #check n_q
  n_q=sapply(n_qs,function(n)
  {
    temp=unique(n)
    if (length(temp)!=1)
    {stop("dimensions mismatch")}
    return((unique(n)))
  })


  #return

  return(list(n_q=n_q,where_q=where))


}

transfo_E <- function(E, type_inter){
  Ecode <- E
  Ecode[which(type_inter=="adj"),2] <- 0
  Ecode[which(type_inter=="diradj"),2] <- -1
  return(Ecode)
}

################################################################################################################
#------------------ transform the data given by the user into  a R6 type coll_interaction object
################################################################################################################

FormattingData = function(listNetwork,vdistrib = NULL)
  #listNet
{

  mats = lapply(listNetwork,function(net){return(net$mat)})
  intFG = lapply(listNetwork,function(net){return(c(net$rowFG,net$colFG))})
  intFG = do.call(rbind,intFG)
  types = sapply(listNetwork,function(net){return(net$type)})


  return(coll_interaction$new(mats,intFG,type = types,distrib = vdistrib))
}


##############################################################################################################
###-----------------------------------FUNCTIONS USED in the VEM algortihm ---------------
##############################################################################################################

#prevent the parameters from being too close from 0 or 1
readjust_alph=function(alpha,eps)
{
  alpha[alpha<eps]=eps
  alpha[alpha>1-eps]=1-eps
  return(alpha/sum(alpha))
}

readjust_pi=function(pi,eps)
{
  pi[pi<eps]=eps
  pi[pi>1-eps]=1-eps
  return(pi)
}


readjust_theta <- function(theta,eps, distrib)
{

  if (distrib == 'bernoulli') {
    theta[theta < eps] = eps
    theta[theta > 1 - eps] = 1 - eps }
  if (distrib == 'poisson') { theta[theta < eps] = eps }

  return(theta)
}


#computing ICL and likelihood
comp_lik_ICL = function(tau,ltheta,lpi,mat_E,list_Mat,n_q,vK,vdistrib)
{
  cardE <- nrow(mat_E)
  Q <-  length(lpi)



  #condlik
  condliks = sapply(1:cardE,function(e)
  {
    gr = mat_E[e,1]
    gc = mat_E[e,2]
    don = list_Mat[[e]]
    if (vdistrib[e] == 'bernoulli') {Unmdon = 1 - don}
    facteur = 1
    if (gc < 1)
    {
      if (gc == 0)  facteur = 1/2 #sbm sym
      gc = gr
      if (vdistrib[e] == 'bernoulli') {diag(Unmdon) = 0}
    }
    if (vdistrib[e] == 'bernoulli') {
      prov = (tau[[gr]]) %*% log(ltheta[[e]]) %*% t(tau[[gc]])
      prov1m = (tau[[gr]]) %*% log(1 - ltheta[[e]]) %*% t(tau[[gc]])
      return((sum(don * prov) + sum(Unmdon * prov1m)) * facteur)
    }
    if (vdistrib[e] == 'poisson') {
      prov = (tau[[gr]]) %*% log(ltheta[[e]]) %*% t(tau[[gc]])
      prov2 = (tau[[gr]]) %*%  ltheta[[e]]  %*% t(tau[[gc]])
      Unit <- matrix(1,nrow = nrow(don),ncol = ncol(don))
      return((sum(don * prov) - sum(Unit * prov2)) * facteur)
    }

  }
  )
  condlik = sum(condliks)

  likmargs = sapply(1:Q,function(q)
  {
    return(sum(tau[[q]]*matrix(log(lpi[[q]]),nrow(tau[[q]]),ncol(tau[[q]]),byrow = TRUE)))
  })
  likmarg = sum(likmargs)

  entros = sapply(1:Q,function(q){return(-sum(log(tau[[q]]) * tau[[q]]))})
  entro = sum(entros)

  #penalty
  pen_mats = vapply(1:cardE,function(s)
  {
    gr = mat_E[s,1]
    gc = mat_E[s,2]
    if (gc > 0) #LBM
    {
      return(c(vK[gr] * vK[gc], prod(dim(list_Mat[[s]]))))
    }
    if (gc == 0) #SBM sym
    {
      return(c(vK[gr] * (vK[gr] + 1) / 2, nrow(list_Mat[[s]]) * (nrow(list_Mat[[s]]) - 1)/2))
    }
    if (gc == -1) #SBM non sym
    {
      return(c(vK[gr] * (vK[gr]), nrow(list_Mat[[s]]) * (nrow(list_Mat[[s]]) - 1)))
    }

  },c(1.0,2.0)
  )


  penICL = sum((vK - 1) * log(n_q))  + sum(pen_mats[1,]) * log(sum(pen_mats[2,]))


  return(list(condlik = condlik,marglik = likmarg,entr = entro,pen = penICL))

}


#distance between pis
distpi=function(pi,pi_old)
{
  Q=length(pi)
  vdis=sapply(1:Q,function(q){
    return(sqrt(sum(as.vector(pi[[q]]-pi_old[[q]])^2)))
  })
  return(sum(vdis))
}

distltheta=function(ltheta,ltheta_old)
{
  Q=length(ltheta)
  vdis=sapply(1:Q,function(q){
    return(sqrt(sum(as.vector(ltheta[[q]]-ltheta_old[[q]])^2)))
  })
  return(sum(vdis))
}




disttau=function(tau,tau_old)
{
  Q=length(tau)
  vdis=sapply(1:Q,function(q){
    return(sqrt(sum(as.vector(tau[[q]]-tau_old[[q]])^2)))
  })
  return(sum(vdis))
}



#verif
ARIS=function(res,classif_vrai)
{
  Q=length(res$tau)
  a=sapply(1:Q,function(q){
    cl=apply(res$tau[[q]],1,which.max)
    adjustedRandIndex(cl,classif_vrai[[q]])})
  return(a)
}

frobenius <- function(X){
  return(sqrt(sum((X)^2)))
}

erparam=function(res,param_vrai)
{
  Q=length(res$alpha)
  cardE=length(res$pi)

  eralpha=sapply(1:Q,function(q){
    C=param_vrai$alpha[[q]]
    alph=res$alpha[[q]]
    K_q=length(res$alpha[[q]])
    o=min(sapply(permn(K_q),
                 function(perm) {
                   A <- alph[perm]
                   return(frobenius(C-A)/frobenius(C))}))
    return(o)
  })

  erpi=sapply(1:cardE,function(q)
  {
    C=param_vrai$pi[[q]]
    pi=res$pi[[q]]
    K_q=nrow(res$pi[[q]])
    K_qprime=ncol(res$pi[[q]])
    o=min(sapply(permn(K_q),
                 function(perm) {
                   oprime=min(sapply(permn(K_qprime),function(perm2){
                     A <- pi[perm,perm2]
                     return(frobenius(C-A)/frobenius(C))}))

                   return(oprime)})
    )

    return(o)

  })

  return(list(eralpha,erpi))
}






##############################################################################################################
# ------------------ Cleaning a list of estimations
##############################################################################################################

Cleaning_estim <- function(dataR6,R){

  #browser()
  ### first  : for equal model keep the estimate wit the highest J
  J_seq <- sapply(R ,function(u){max(u$vJ)})
  o <- order(J_seq,decreasing = TRUE)
  R.ordered.1 <- lapply(o,function(i){R[[i]]})

  if (dataR6$Q == 1) {
  seq_nb_clust <- cbind((sapply(R.ordered.1,function(u){u$param_estim$vK})),J_seq[o],1:length(R))
  } else {
  seq_nb_clust <- cbind(t(sapply(R.ordered.1,function(u){u$param_estim$vK})),J_seq[o],1:length(R))
  }


  if (length(R) > 1)
  {
    seq_nb_clust <- seq_nb_clust[!duplicated(seq_nb_clust[,1:dataR6$Q]),]
    if (is.null(nrow(seq_nb_clust))) {seq_nb_clust = matrix(seq_nb_clust,nrow=1)}
    R <-  R.ordered.1[seq_nb_clust[,dataR6$Q + 2]]
  }

  ### Order the results by ICL
  ICL_seq <- sapply(R,function(u){u$ICL})
  o <- order(ICL_seq,decreasing = TRUE)
  res <- lapply(o,function(i){R[[i]]})
  return(res)
}


compar_classif <- function(classif_1,classif_2){

  if (length(classif_1) != length(classif_2) )
  {
    stop('classification can not be compared because objects of different sizes')
  }
  ARI <- mean(vapply(1:length(classif_1),function(q){adjustedRandIndex(classif_1[[q]],classif_2[[q]])},1))
  if (ARI == 1){return(TRUE)}else{return(FALSE)}
}

clean_collection_classif <- function(collection_classif,ind_ref){

  L <- length(collection_classif)
  mat_compar <- matrix(NA,L,L)
  for (l in 2:L){for (k in 1:(l - 1)){mat_compar[l,k] <- compar_classif(collection_classif[[k]],collection_classif[[l]])}}
  w <- which(rowSums(mat_compar[ind_ref:L,],na.rm = TRUE)>0)- ind_ref+1
  u <- (1:L)[ -w]
  res <- lapply(u,function(i){collection_classif[[i]]})
  return(res)
}




  }





}




