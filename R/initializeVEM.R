initialize = function(dataR6 , param , method = "CAH", classifInit = NULL)
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
    groups <- lapply(1:Q,function(q){cahClassif(dataR6, q, v_K[q])})
    }

  if (method == "given") { groups <- classifInit}


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
      gr1 <- cahClassif(dataR6, q, 2, subset = which(classif_q == k)); gr <- gr1; gr[gr1 == 1] <- k; gr[gr1 == 2] <- K_q + 1
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







#---------------------------------------------------------------------------------
# NA-safe distances for the hierarchical clustering
#---------------------------------------------------------------------------------

# Manhattan distance between the rows (side = 1) or the columns (side = 2) of mat.
# stats::dist already ignores the NA and rescales, but returns NA when two
# individuals have no observed entry in common: these are replaced by the mean distance.
distNA <- function(mat, side = 1, subset = NULL) {
  if (side == 2) {mat <- t(mat)}
  if (!is.null(subset)) {mat <- mat[subset, , drop = FALSE]}
  d <- stats::dist(mat, "manhattan")
  if (anyNA(d)) {
    m <- mean(d, na.rm = TRUE)
    d[is.na(d)] <- if (is.finite(m)) m else 0
  }
  return(d)
}

# Distance between the individuals of FG q, summed over all the networks involving q
distFG <- function(dataR6, q, subset = NULL) {
  w_q <- dataR6$where[[q]]
  dists <- lapply(seq_len(nrow(w_q)), function(i) {
    distNA(dataR6$mats[[w_q[i, 1]]], side = w_q[i, 2], subset = subset)
  })
  return(Reduce('+', dists))
}

# Ward hierarchical clustering of the individuals of FG q (or of a subset of them) in K blocks
cahClassif <- function(dataR6, q, K, subset = NULL) {
  n <- if (is.null(subset)) dataR6$v_NQ[q] else length(subset)
  if (K <= 1 || n <= 1) {return(rep(1, n))}
  cah <- stats::hclust(distFG(dataR6, q, subset), method = "ward.D")
  return(stats::cutree(cah, min(K, n)))
}

# Modify the classification classif_q of FG q so that it has exactly K blocks:
# the biggest block is split in two (CAH) or the two closest blocks (average distance) are merged.
adjustClassif <- function(classif_q, K, dataR6, q) {
  classif_q <- match(classif_q, unique(sort(classif_q)))
  K_q <- max(classif_q)
  if (K_q == K) {return(classif_q)}
  D <- as.matrix(distFG(dataR6, q))
  while (K_q > K) {
    couples <- utils::combn(1:K_q, 2)
    dmean <- apply(couples, 2, function(kk) {
      mean(D[classif_q == kk[1], classif_q == kk[2]])
    })
    best <- couples[, which.min(dmean)]
    classif_q[classif_q == best[2]] <- best[1]
    classif_q <- match(classif_q, unique(sort(classif_q)))
    K_q <- K_q - 1
  }
  while (K_q < K) {
    nk <- tabulate(classif_q, K_q)
    if (max(nk) < 2) {break}
    k <- which.max(nk)
    ind <- which(classif_q == k)
    gr <- cahClassif(dataR6, q, 2, subset = ind)
    classif_q[ind[gr == 2]] <- K_q + 1
    K_q <- K_q + 1
  }
  return(classif_q)
}


#---------------------------------------------------------------------------------
# Initialisation with blockmodels (SBM or LBM fitted on each network separately)
#---------------------------------------------------------------------------------

# Fit a SBM / LBM on network e with blockmodels. Returns a list with the FG indices
# and the corresponding classifications, or NULL if the distribution is not handled.
blockmodelsOneNet <- function(dataR6, e, v_Kmin, v_Kmax, verbose = FALSE, nbCores = NULL) {
  bm_model <- switch(dataR6$v_distrib[e],
    "bernoulli" = blockmodels::BM_bernoulli,
    "poisson" = blockmodels::BM_poisson,
    "gaussian" = blockmodels::BM_gaussian,
    NULL
  )
  if (is.null(bm_model)) {
    if (verbose) {
      print(paste0("blockmodels initialisation not available for the ", dataR6$v_distrib[e],
                   " distribution: network ", e, " is not used for the initialisation"))
    }
    return(NULL)
  }
  isInc <- dataR6$typeInter[e] == "inc"
  indFG <- if (isInc) dataR6$E[e, ] else dataR6$E[e, 1]
  membership_type <- if (isInc) "LBM" else if (dataR6$typeInter[e] == "diradj") "SBM" else "SBM_sym"
  ncores <- if (is.null(nbCores) || .Platform$OS.type == "windows") 1 else nbCores

  estim <- bm_model(
    membership_type = membership_type,
    adj = dataR6$mats[[e]],
    verbosity = 0,
    plotting = "",
    explore_max = sum(v_Kmax[indFG]),
    ncores = ncores
  )
  estim$estimate()

  # classification(s) for each explored model
  getZ <- function(memb) {
    Zs <- if (isInc) list(memb$Z1, memb$Z2) else list(memb$Z)
    lapply(Zs, function(Z) {z <- max.col(Z, ties.method = "first"); match(z, unique(sort(z)))})
  }
  ICL <- estim$ICL
  allZ <- lapply(seq_along(ICL), function(i) {
    if (is.null(estim$memberships[[i]]) || !is.finite(ICL[i])) NULL else getZ(estim$memberships[[i]])
  })
  # keep the best ICL among models compatible with v_Kmin, v_Kmax
  admissible <- vapply(allZ, function(Z) {
    !is.null(Z) && all(vapply(Z, max, 1) >= v_Kmin[indFG]) && all(vapply(Z, max, 1) <= v_Kmax[indFG])
  }, TRUE)
  ICLadm <- ifelse(admissible, ICL, -Inf)
  if (all(!is.finite(ICLadm))) {ICLadm <- ifelse(vapply(allZ, is.null, TRUE), -Inf, ICL)}
  if (all(!is.finite(ICLadm))) {return(NULL)}
  return(list(indFG = indFG, Z = allZ[[which.max(ICLadm)]]))
}


# Build the collection of initial classifications from blockmodels fits on each network.
# For each FG, every network involving it gives a candidate classification; all the combinations are returned.
# If v_K is given, each candidate is adjusted to have exactly v_K[q] blocks; otherwise
# candidates are only forced within [v_Kmin, v_Kmax].
# FGs not covered by any network handled by blockmodels are initialised by CAH.
initBlockmodels <- function(dataR6, v_Kmin, v_Kmax, v_K = NULL, verbose = FALSE, nbCores = NULL) {
  Q <- dataR6$Q
  candidates <- lapply(1:Q, function(q) {list()})
  for (e in 1:dataR6$cardE) {
    res_e <- blockmodelsOneNet(dataR6, e, v_Kmin, v_Kmax, verbose = verbose, nbCores = nbCores)
    if (!is.null(res_e)) {
      for (j in seq_along(res_e$indFG)) {
        q <- res_e$indFG[j]
        candidates[[q]] <- c(candidates[[q]], list(res_e$Z[[j]]))
      }
    }
  }

  for (q in 1:Q) {
    if (length(candidates[[q]]) == 0) {
      Kq <- if (is.null(v_K)) v_Kmin[q] else v_K[q]
      candidates[[q]] <- list(cahClassif(dataR6, q, Kq))
    }
    candidates[[q]] <- lapply(candidates[[q]], function(z) {
      Kz <- max(z)
      Ktarget <- if (!is.null(v_K)) v_K[q] else min(max(Kz, v_Kmin[q]), v_Kmax[q])
      z <- adjustClassif(z, Ktarget, dataR6, q)
      names(z) <- dataR6$namesInd[[q]]
      z
    })
    # remove duplicated candidates (equal up to label switching)
    keep <- rep(TRUE, length(candidates[[q]]))
    for (i in seq_along(candidates[[q]])[-1]) {
      for (j in 1:(i - 1)) {
        if (keep[j] && samePartition(candidates[[q]][[i]], candidates[[q]][[j]])) {keep[i] <- FALSE}
      }
    }
    candidates[[q]] <- candidates[[q]][keep]
  }

  combin <- as.matrix(expand.grid(lapply(candidates, seq_along)))
  collection <- lapply(seq_len(nrow(combin)), function(i) {
    lapply(1:Q, function(q) candidates[[q]][[combin[i, q]]])
  })
  return(collection)
}

# TRUE if the two classifications are equal up to label switching
samePartition <- function(z1, z2) {
  length(unique(z1)) == length(unique(z2)) &&
    length(unique(paste(z1, z2))) == length(unique(z1))
}
