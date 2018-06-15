initialize=function(data,param,method="CAH",givenclassif=NULL)
  #method can be pure random or classical clustering CAH or given classif
{
  mats <- data$mats;
  Ecode <- data$Ecode; 
  vK <- param$vK;
  v_NQ <- data$v_NQ;
  where_q <- data$where;
  Q <- data$Q
  card_E <- data$card_E
  #parameters
  eps=.Machine$double.eps
  
  if (method=="random"){groups=lapply(1:Q,function(q){return(sample(1:vK[q],v_NQ[q],replace=TRUE))})}
  
  if (method=="CAH")
  {
    #compute distances between individuals
    groups=lapply(1:Q,function(q)
    {
      w_q=where_q[[q]]
      dists=lapply(as.list(as.data.frame(t(w_q))),function(l) 
      { 
        if (l[2]==1)   d=dist(mats[[l[1]]],"manhattan")
        else  d=dist(t(mats[[l[1]]]),"manhattan")
        return(d)
      })
      totdist=Reduce('+',dists)
      
      #if more than 1 element
       if (vK[q]>1)
       {
        cah1=hclust(totdist,method="ward.D")
        gr1=cutree(cah1,vK[q])  
       }
       else {gr1=rep(1,data$v_NQ[q])}
      
      return(gr1)
    })
    }
  
  if (method=="given"){groups=givenclassif}
  
  tau=  lapply(1:Q,function(j){
    matc=matrix(eps,v_NQ[j],vK[j],byrow=TRUE)
    matc[cbind(1:v_NQ[j],groups[[j]])]=1-eps
    #normalize tau
    matc=matc/rowSums(matc)
    return(matc)
  })
  return(list(groups=groups,tau=tau))
}  


###################################################################"
sequential_initialize=function(classif.c, data ,Kmin=NULL,Kmax=NULL,os = "windows"){

  Q <- length(classif.c); 
  if(is.null(Kmin)){Kmin  = rep(1,Q)}
  if(is.null(Kmax)){Kmax  = rep(10,Q)}
  if(length(Kmin)==1){Kmin = rep(Kmin,Q)}
  if(length(Kmax)==1){Kmax = rep(Kmax,Q)}
  
  F_q <-  function(q){
    classif_forward_q  <- classif_backward_q<-  classif.c;
    classif_forward_q <- split_classif(classif.c,q,data,Kmax[q])
    classif_backward_q <- merge_classif(classif.c,q,Kmin[q])
    res <- do.call(c, list(classif_backward_q, classif_forward_q))
    return(res)}
  
  list_classif_init = list()
  for (q in 1:Q){list_classif_init <- do.call(c,list(list_classif_init,F_q(q)))}
  return(list_classif_init)
}
# 
vK_init=c(2,2,2)

calc_vK = function(classif){
  Q = length(classif); 
  vK <- sapply(1:Q,function(q){length(unique(classif[[q]]))})
  return(vK)
}
  


split_classif =  function(classif,q,data,Kmax_q){
  
  classif_q <- classif[[q]]
  classif_q = match(classif_q, unique(sort(classif_q))) # pour avoir comme numéro de groupe 1... K_q
  K_q <- max(classif_q)
  if(K_q < Kmax_q ){
    nq = table(classif_q)
    uq = unique(sort(classif_q))
    wq = uq[which(nq>1)]
    new_classif <-  lapply(wq,function(k){
      classif_split_q_k <-  classif; # nouvelle classif coupant le cluster k du FG q en 2
      
      # on coupe les individus du groupe $k$ en utilisant CAH sur les distances
      dists=lapply(as.list(as.data.frame(t(data$where[[q]]))),function(l) 
        { if (l[2]==1)   d=dist(data$mats[[l[1]]][classif_q==k,],"manhattan")
          else  d=dist(t(data$mats[[l[1]]][,classif_q==k]),"manhattan")
          return(d)
        })
      totdist=Reduce('+',dists)
      cah1=hclust(totdist,method="ward.D")
      gr1=cutree(cah1,2); gr <-gr1; gr[gr1==1] = k; gr[gr1==2] = K_q+1
      classif_split_q_k[[q]][classif_q==k] <- gr #nqk = sum(classif_q==k); sample(c(k,K_q+1),nqk,replace=TRUE)
      
      if(length(unique(gr))==1){
        stop("parfois l'initilisation ")
      }
      
      
      return(classif_split_q_k)})
  }
  else{new_classif=list()}
  return(new_classif)
}



merge_classif =  function(classif,q,Kmin_q){
  
  classif_q <- classif[[q]]
  classif_q = match(classif_q, unique(sort(classif_q)))
  K_q <- max(classif_q)
  if(K_q > Kmin_q){
    couples_q <- t(combn(1:K_q, 2))
    R = K_q*(K_q-1)/2
    new_classif<-  lapply(1:R,function(k){
      classif_merge_q_k1k2 = classif;   # nouvelle classif réunissant les cluster k1 et k2 du FG q
      k1 <- couples_q[k,1]
      k2 <- couples_q[k,2]
      z_q_k1k2 <- classif_q;
      z_q_k1k2[classif_q==k2] <- k1;
      # renumbering of the classification to get 1:K_q-1
      classif_merge_q_k1k2[[q]] = match(z_q_k1k2, unique(sort(z_q_k1k2)))
      return(classif_merge_q_k1k2)})
  }
  else{new_classif=list()}
  return(new_classif)
}
#########################################################################

  
  
  
  
  
  