
#initialization provides a first classification and corresponding tau
#method can be random, CAH, given in which case add given classif

varEMMBM <- function(dataR6,classifInit,tauInit=NULL)
  #data :  coll_interaction type object
  #classif  : liste de classifications (Z) au sein des groupes fonctionnels fg (de longueur Q )

{
  #checking the dimensions of matrices and extracting number of indiviuals
  #for a given q functional groups, gives the list of matrix in row or in columns where it plays a role


  where_q <- dataR6$where
  n_q <- dataR6$v_NQ
  cardE <- dataR6$cardE
  matE <- dataR6$Ecode ### a remettre en public?????
  list_Mat <- dataR6$mats
  v_distrib <- dataR6$v_distrib

  #const_lik_poisson <- vapply(1:cardE, function(e){if (v_distrib[e] == 'poisson') {return(sum(lfactorial(c(list_Mat[[e]]))))} else {return(0)}},1)
  #const_lik_poisson <- sum(const_lik_poisson)

  #some constants
  eps <- 2*.Machine$double.eps
  valStopCrit <- 1e-6



  ##  initialisation
  v_K <- calcVK(classifInit)
  tau <- tauInit
  if (is.null(tau)) {
    tau <-  lapply(1:dataR6$Q,function(j){
      mat <- matrix(eps,n_q[j],v_K[j],byrow = TRUE)
      mat[cbind(1:n_q[j],classifInit[[j]])] <- 1 - eps
      #normalize tau
      mat <- mat/rowSums(mat)
      return(mat)
    })
  }



  #for computations when SBM putting 0 on the diagonal
  indSBM <- which(matE[,2] < 1)
  for (i in indSBM) {diag(list_Mat[[i]])  <- 0}


  #entering VEM
  maxiter <- 1000
  stopcrit <- 0
  iterVEM <- 0

  vJ  <- numeric(maxiter)

  #for stopping criterion
  list_theta <- lapply(1:cardE,function(s)
  {
    gr <- matE[s,1]
    gc <- matE[s,2]
    if (gc < 1) gc <- gr
    return(matrix(Inf,v_K[gr],v_K[gc]))
  })




  ######################
  # Algo begins
  #####################
  while (iterVEM < maxiter & stopcrit == 0)
  {
    iterVEM <- iterVEM + 1


    #--------------------------------   M step
    list_theta_old <- list_theta

    list_pi = lapply(tau,colMeans)
    list_theta  = lapply(1:cardE,function(j){
      gr <- matE[j,1]
      gc <- matE[j,2]


      if (gc < 1) {  #for sbm sym or notsym
        gc <- gr
        #useful matrix
        Unitmdiag <- matrix(1,nrow = n_q[gr],ncol = n_q[gc])
        diag(Unitmdiag) <- 0
        #bernoulli or poisson distribution same expression for M step
        #list_thetac = t(tau[[gr]]) %*% list_Mat[[j]] %*% tau[[gc]] / (t(tau[[gr]]) %*% (Unitmdiag) %*% tau[[gc]])
        #   passage en crossprod
        list_thetac <- crossprod(crossprod(list_Mat[[j]], tau[[gr]]), tau[[gc]]) / crossprod(crossprod(Unitmdiag, tau[[gr]]), tau[[gc]])
      }
      else #for lbm
      {
        Unit <- matrix(1,nrow = n_q[gr],ncol = n_q[gc])
        #list_thetac <-  t(tau[[gr]]) %*% list_Mat[[j]] %*% tau[[gc]] / (t(tau[[gr]]) %*% (Unit) %*% tau[[gc]])
        list_thetac <- crossprod(crossprod(list_Mat[[j]], tau[[gr]]), tau[[gc]]) / crossprod(crossprod(Unit, tau[[gr]]), tau[[gc]])

      }
      return(list_thetac)})

    #prevent the values from being too close from 0 or 1
    list_pi <- lapply(list_pi,readjustPi,eps)
    list_theta <- lapply(1:cardE,function(e){readjustTheta(list_theta[[e]],eps,v_distrib[e])})


    # pseudolik <- comp_lik_ICL(tau,list_theta,list_pi,matE,list_Mat,n_q,v_K)
    # vJ <- pseudolik$condlik + pseudolik$marglik + pseudolik$entr
    # print(vJ)

    #stop criterion
    if (distListTheta(list_theta,list_theta_old) < valStopCrit) stopcrit <- 1


    #--------------------------------   VE step : boucle
    iterVE <- 0
    stopVE <- 0
    while ((iterVE < maxiter) & (stopVE == 0))
    {
      #VE step

      tau_old <- tau #useful ?
      #browser()
      for (q in 1:dataR6$Q)
      {


        if (v_K[q] == 1) { tau[[q]] = matrix(1,ncol  = 1,nrow = n_q[q])}
        else{
          w_q <- where_q[[q]]
          der <- lapply(as.list(as.data.frame(t(w_q))),function(l)
            {
            second_index <- matE[l[1],ifelse(l[2] == 1,2,1)]
            qprime <- second_index

            don <- list_Mat[[l[1]]] #pb a regler ligne ou colonne et si sbm enlever de 1-don la diag
            if (v_distrib[l[1]] == 'bernoulli' ) { Unmdon <- 1 - don }

            matlist_theta <- list_theta[[l[1]]]

            if (qprime < 1)   #if sbm
            {
              if (v_distrib[l[1]] == 'bernoulli' ) {diag(Unmdon) <- 0 }
              qprime <- q
            }
            else  # if lbm
            {
              if (l[2] == 2) #functional group q at stake in rows
              {
                don <- t(don)
                if (v_distrib[l[1]] == 'bernoulli' ) {Unmdon <- t(Unmdon) }
                matlist_theta <- t(matlist_theta)
              }
            }
            #browser()
            #poisson or bernoulli likelihood
            switch(v_distrib[l[1]],
              bernoulli = {lik  = don %*% tcrossprod(tau[[qprime]],log(matlist_theta)) + Unmdon %*% tcrossprod(tau[[qprime]],log(1 - matlist_theta))},
              poisson = {Unit <- matrix(1,nrow(don),ncol(don))
              lik  <- -Unit %*% tcrossprod(tau[[qprime]], matlist_theta) + don %*% tcrossprod(tau[[qprime]],log(matlist_theta))  }
            )



            if (second_index < 0)#sbm but non sym
            {
              don <- t(don)
              if (v_distrib[l[1]] == 'bernoulli' ) { Unmdon <- t(Unmdon) }
              matlist_theta <- t(matlist_theta)
              switch(v_distrib[l[1]],
                #bernoulli = {lik = lik + don %*% tau[[qprime]] %*% t(log(matlist_theta)) + Unmdon %*% tau[[qprime]] %*% t(log(1 - matlist_theta))},
                bernoulli = {lik = lik + don %*% tcrossprod(tau[[qprime]],log(matlist_theta)) + Unmdon %*% tcrossprod(tau[[qprime]],log(1 - matlist_theta))},
                poisson = {lik = lik + don %*% tcrossprod(tau[[qprime]],log(matlist_theta))  -  matrix(1,nrow(don),ncol(don)) %*% tcrossprod(tau[[qprime]], matlist_theta)})
            }
            return(lik)
          })# ----- fin der
          L <- (Reduce('+',der))

          B <- L + matrix(log(list_pi[[q]]),nrow = nrow(tau[[q]]),ncol = v_K[q],byrow = TRUE)
          #B <- B - max(B)
          B <- B - matrix(apply(B,1,mean),nrow = n_q[q],ncol = v_K[q],byrow = FALSE)
          B[B > 709] = 709


          temp <- exp(B)
          temp2 <- temp/rowSums(temp)
          temp2[temp2 < eps] <- eps
          temp2[temp2 > (1 - eps)] <- 1 - eps
          temp2 <- temp2/rowSums(temp2)
          #if (any(is.na(temp2))) { browser()}
          #if (any(is.infinite(temp2))) { browser()}

          tau[[q]] <- temp2


        }
        #if (any(dim(tau[[q]]) != dim(tau_old[[q]]))) { browser()}
      }
      iterVE <- iterVE + 1
      if (distTau(tau,tau_old) < valStopCrit)   stopVE <- 1
    }#-------------------- END of VE Step


    #computing lik
    pseudolik <- compLikICL(tau,list_theta,list_pi,matE,list_Mat,n_q,v_K,v_distrib)
    vJ[iterVEM] <- pseudolik$condLik + pseudolik$margLik + pseudolik$entr
  } # ------------ end of EM var

  #computing ICL
  likICL <- compLikICL(tau,list_theta,list_pi,matE,list_Mat,n_q,v_K,v_distrib)

 ICL <-  likICL$condLik + likICL$margLik - 1/2 * likICL$pen

  paramEstim   <- MBMfit$new(v_K = v_K, v_distrib = v_distrib, list_pi = list_pi,list_theta = list_theta);
  paramEstim$tau <- tau
  vJ <- vJ[1:iterVEM]
  return(list(paramEstim = paramEstim,ICL = ICL,vJ = vJ))
}

