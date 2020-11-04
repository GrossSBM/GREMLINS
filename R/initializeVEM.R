initialize = function(dataR6 , param , method = "CAH", givenclassif = NULL)
  #method can be pure random or classical clustering CAH or given classif
{

  v_K <- param$v_K;
  v_NQ <- dataR6$v_NQ;
  where_q <- dataR6$where;
  Q <- dataR6$Q
  eps = .Machine$double.eps

  if (method == "random") { groups <- lapply(1:Q,function(q){return(sample(1:v_K[q],v_NQ[q],replace = TRUE))})}

  if (method == "CAH")
  {
    #compute distances between individuals
    groups <- lapply(1:Q,function(q)
    {
      w_q <- where_q[[q]]
      dists <- lapply(as.list(as.data.frame(t(w_q))),function(l)
      {
        if (l[2] == 1)   d <- stats::dist(dataR6$mats[[l[1]]],"manhattan")
        else  d <- stats::dist(t(dataR6$mats[[l[1]]]),"manhattan")
        return(d)
      })
      totdist <- Reduce('+',dists)

      #if more than 1 element
       if (v_K[q] > 1)
       {
        cah1 <- stats::hclust(totdist,method = "ward.D")
        gr1 <- stats::cutree(cah1,v_K[q])
       }
       else {gr1 <- rep(1,dataR6$v_NQ[q])}

      return(gr1)
    })
    }

  if (method == "given") { groups <- givenclassif}


  tau <-  lapply(1:Q,function(j){
    matc <- matrix(eps,v_NQ[j],v_K[j],byrow = TRUE)
    matc[cbind(1:v_NQ[j],groups[[j]])] <- 1 -  eps
    #normalize tau
    matc <- matc/rowSums(matc)
    return(matc)
  })
  return(list(groups = groups,tau = tau))
}


#---------------------------------------------------------------------------------
sequentialInitialize = function(classif.c, dataR6 ,Kmin = NULL, Kmax = NULL, os = "windows"){

  Q <- length(classif.c);
  if (is.null(Kmin)) {Kmin  <- rep(1,Q)}
  if (is.null(Kmax)) {Kmax  <- rep(10,Q)}
  if (length(Kmin) == 1) {Kmin <- rep(Kmin,Q)}
  if (length(Kmax) == 1) {Kmax <- rep(Kmax,Q)}

  F_q <-  function(q){
    classifForward_q  <- classifBackward_q <-  classif.c;
    classifForward_q <- splitClassif(classif.c,q,dataR6,Kmax[q])
    classifBackward_q <- mergeClassif(classif.c,q,Kmin[q])
    res <- do.call(c, list(classifBackward_q, classifForward_q))
    return(res)}

  listClassifInit <- list()
  for (q in 1:Q) { listClassifInit <- do.call(c,list(listClassifInit,F_q(q)))}
  return(listClassifInit)
}
#




#---------------------------------------------------------------------------------

calcVK = function(classif){
  Q = length(classif);
  v_K <- sapply(1:Q,function(q){length(unique(classif[[q]]))})
  return(v_K)
}


#---------------------------------------------------------------------------------

splitClassif =  function(classif,q,dataR6,Kmax_q){

  classif_q <- classif[[q]]
  classif_q <- match(classif_q, unique(sort(classif_q))) # pour avoir comme numéro de groupe 1... K_q
  K_q <- max(classif_q)
  if (K_q < Kmax_q ) {
    nq <- table(classif_q)
    uq <- unique(sort(classif_q))
    wq <- uq[which(nq > 1)]
    newClassif <-  lapply(wq,function(k){
      classif_split_q_k <-  classif; # nouvelle classif coupant le cluster k du FG q en 2

      # on coupe les individus du groupe $k$ en utilisant CAH sur les distances
      dists <- lapply(as.list(as.data.frame(t(dataR6$where[[q]]))),function(l)
        { if (l[2] == 1)   d <- stats::dist(dataR6$mats[[l[1]]][classif_q == k,],"manhattan")
          else  d <- stats::dist(t(dataR6$mats[[l[1]]][,classif_q == k]),"manhattan")
          return(d)
        })
      totdist <- Reduce('+',dists)
      cah1 <- stats::hclust(totdist,method = "ward.D")
      gr1 <- stats::cutree(cah1,2); gr <- gr1; gr[gr1 == 1] <- k; gr[gr1 == 2] <- K_q + 1
      classif_split_q_k[[q]][classif_q == k] <- gr #nqk = sum(classif_q==k); sample(c(k,K_q+1),nqk,replace=TRUE)

      return(classif_split_q_k)})
  }
  else{newClassif <- list()}
  return(newClassif)
}


#---------------------------------------------------------------------------------

mergeClassif =  function(classif,q,Kmin_q){

  classif_q <- classif[[q]]
  classif_q = match(classif_q, unique(sort(classif_q)))
  K_q <- max(classif_q)
  if (K_q > Kmin_q) {
    couples_q <- t(utils::combn(1:K_q, 2))
    R = K_q * (K_q - 1)/2
    newClassif <-  lapply(1:R,function(k){
      classif_merge_q_k1k2 = classif;   # nouvelle classif réunissant les cluster k1 et k2 du FG q
      k1 <- couples_q[k,1]
      k2 <- couples_q[k,2]
      z_q_k1k2 <- classif_q;
      z_q_k1k2[classif_q == k2] <- k1;
      # renumbering of the classification to get 1:K_q-1
      classif_merge_q_k1k2[[q]] = match(z_q_k1k2, unique(sort(z_q_k1k2)))
      return(classif_merge_q_k1k2)})
  }
  else{newClassif = list()}
  return(newClassif)
}

#---------------------------------------------------------------------------------





