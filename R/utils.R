###------------------------------------------------------------#####

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




###-----------------------------------FUNCTIONS USED in the VEM algortihm ---------------

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


readjust_theta=function(theta,eps)
{
  theta[theta<eps]=eps
  theta[theta>1-eps]=1-eps
  return(theta)
}


#computing ICL and likelihood
comp_lik_ICL=function(tau,pi,alpha,mat_E,list_Mat,n_q,vK)
{
  cardE=nrow(mat_E)
  Q=length(alpha)

  #condlik
  condliks=sapply(1:cardE,function(e)
  {
    gr=mat_E[e,1]
    gc=mat_E[e,2]
    don=list_Mat[[e]]
    Unmdon=1-don
    facteur=1
    if (gc<1)
    {
      if (gc==0)  facteur=1/2 #sbm sym
      gc=gr
      diag(Unmdon)=0
    }
    prov=(tau[[gr]])%*%log(pi[[e]])%*%t(tau[[gc]])
    prov1m=(tau[[gr]])%*%log(1-pi[[e]])%*%t(tau[[gc]])
    return((sum(don*prov) +sum(Unmdon*prov1m))*facteur)
  }
  )
  condlik=sum(condliks)

  likmargs=sapply(1:Q,function(q)
  {
    return(sum(tau[[q]]*matrix(log(alpha[[q]]),nrow(tau[[q]]),ncol(tau[[q]]),byrow = TRUE)))
  })
  likmarg=sum(likmargs)

  entros=sapply(1:Q,function(q)
  {
    return(-sum(log(tau[[q]])*tau[[q]]))
  })

  entro = sum(entros)

  #penalty
  pen_mats=vapply(1:cardE,function(s)
  {
    gr=mat_E[s,1]
    gc=mat_E[s,2]
    if (gc>0) #LBM
    {
      return(c(vK[gr]*vK[gc], prod(dim(list_Mat[[s]]))))
    }
    if (gc==0) #SBM sym
    {
      return(c(vK[gr]*(vK[gr]+1)/2, nrow(list_Mat[[s]])*(nrow(list_Mat[[s]])-1)/2))
    }
    if (gc==-1) #SBM non sym
    {
      return(c(vK[gr]*(vK[gr]), nrow(list_Mat[[s]])*(nrow(list_Mat[[s]])-1)))
    }

  },c(1.0,2.0)
  )


  penICL=sum((vK-1)*log(n_q))  + sum(pen_mats[1,])*log(sum(pen_mats[2,]))


  return(list(condlik=condlik,marglik=likmarg,entr=entro,pen=penICL))

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





# internal function -------------------------------------------------------

FormattingData = function(listNetwork)
#listNet
{
  nNet = length(listNetwork)
  mats = lapply(listNetwork,function(net){return(net$mat)})
  intFG = lapply(listNetwork,function(net){return(c(net$rowFG,net$colFG))})
  intFG = do.call(rbind,intFG)
  types = sapply(listNetwork,function(net){return(net$type)})

  return(coll_interaction$new(mats,intFG,type = types))
}



