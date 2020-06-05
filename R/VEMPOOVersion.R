
#initialization provides a first classification and corresponding tau
#method can be random, CAH, given in which case add given classif

varEMMBM <- function(dataR6,classifInit,tauInit = NULL, maxiterVE = NULL,maxiterVEM = NULL)
  #data :  coll_interaction type object
  #classif  : liste de classifications (Z) au sein des groupes fonctionnels fg (de longueur Q )

{
  #checking the dimensions of matrices and extracting number of indiviuals
  #for a given q functional groups, gives the list of matrix in row or in columns where it plays a role



  namesFG <- dataR6$namesFG
  where_q <- dataR6$where
  n_q <- dataR6$v_NQ
  cardE <- dataR6$cardE
  matE <- dataR6$Ecode ### a remettre en public?????
  #list_Mat <- dataR6$mats
  list_MaskNA <- lapply(dataR6$mats,function(m){1 * (1 - is.na(m))})
  list_Mat <- lapply(dataR6$mats,function(m){m[is.na(m)] = 0; m})

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
  #indSBM <- which(matE[,2] < 1)
  #for (i in indSBM) {diag(list_Mat[[i]])  <- 0}


  #entering VEM
  if (is.null(maxiterVE)) { maxiterVE = 1000}
  if (is.null(maxiterVEM)) { maxiterVEM  =  1000}
  maxiter <- maxiterVEM
  stopcrit <- 0
  iterVEM <- 0

  vJ  <- numeric(maxiter)

  #for stopping criterion
  list_theta <- lapply(1:cardE,function(e)
  {
    gr <- matE[e,1]
    gc <- matE[e,2]
    if (gc < 1) gc <- gr
    if (v_distrib[e] == 'gaussian') {return(list(mean = matrix(Inf,v_K[gr],v_K[gc]), var = matrix(Inf,v_K[gr],v_K[gc])))}
    if (v_distrib[e] == 'ZIgaussian') {
      return(list(mean = matrix(Inf,v_K[gr],v_K[gc]), var = matrix(Inf,v_K[gr],v_K[gc]),p0 = matrix(Inf,v_K[gr],v_K[gc])))}
    if (v_distrib[e] %in% c('bernoulli','poisson','laplace')) { return(matrix(Inf,v_K[gr],v_K[gc]))}
  }
  )




  ######################
  # Algo begins
  #####################
  no.convergence = 0;
  while (iterVEM < maxiter & stopcrit == 0)
  {
    iterVEM <- iterVEM + 1


    #--------------------------------   M step
    list_theta_old <- list_theta

    list_pi = lapply(tau,colMeans)
    list_theta  = lapply(1:cardE,function(e){
      gr <- matE[e,1]
      gc <- matE[e,2]
      if (gc < 1) {  #for sbm sym or notsym
         gc <- gr
      }
      #   #useful matrix
      #   Unitmdiag <- matrix(1,nrow = n_q[gr],ncol = n_q[gc])
      #   diag(Unitmdiag) <- 0
      #   Unit <- Unitmdiag
      # }else{
      #   Unit <- matrix(1,nrow = n_q[gr],ncol = n_q[gc])
      # }

      Denom <- crossprod(crossprod(list_MaskNA[[e]], tau[[gr]]), tau[[gc]])
      mu <- crossprod(crossprod(list_Mat[[e]], tau[[gr]]), tau[[gc]])

      if (v_distrib[e] %in% c('poisson','bernoulli')) { #bernoulli or poisson distribution same expression for M step
        list_theta_e <- mu  / Denom
      }
      if (v_distrib[e]  %in% c('gaussian', 'ZIgaussian')) {
        list_theta_e <- list()
        list_theta_e$mean <- mu / Denom
        A <- crossprod(crossprod(list_Mat[[e]]^2, tau[[gr]]), tau[[gc]]) / Denom
        list_theta_e$var <-  A - list_theta_e$mean^2
      }
      if (v_distrib[e] == 'ZIgaussian') {
        NonZeros_e  <- list_Mat[[e]] != 0
        list_theta_e$mean <- mu / crossprod(crossprod(NonZeros_e, tau[[gr]]), tau[[gc]])
        A <- crossprod(crossprod(list_Mat[[e]]^2, tau[[gr]]), tau[[gc]]) /  crossprod(crossprod(NonZeros_e, tau[[gr]]), tau[[gc]])
        list_theta_e$var <-  A - list_theta_e$mean^2
        list_theta_e$p0 <- 1 - crossprod(crossprod(NonZeros_e, tau[[gr]]), tau[[gc]]) / Denom
      }
      if (v_distrib[e] == 'laplace') {
          Omega <- list_Mat[[e]]
          diag(Omega) <- 0
          list_theta_e <- crossprod(crossprod(abs(Omega), tau[[gr]]), tau[[gc]]) / Denom
      }
      return(list_theta_e)})

    #prevent the values from being too close from 0 or 1 for pi, and from the bounds of the support for the others
    list_pi <- lapply(list_pi,readjustPi,eps)
    list_theta <- lapply(1:cardE,function(e){readjustTheta(list_theta[[e]],eps,v_distrib[e])})



    #stop criterion
    if (distListTheta(list_theta,list_theta_old) < valStopCrit) stopcrit <- 1


    #--------------------------------   VE step : boucle
    iterVE <- 0
    stopVE <- 0

    while ((iterVE < maxiterVE) & (stopVE == 0))
    {
      #VE step

      tau_old <- tau #useful ?

        for (q in 1:dataR6$Q)
          {
          #----------------------------------------------------------------
          if (v_K[q] == 1) { tau[[q]] = matrix(1,ncol  = 1,nrow = n_q[q])
          }else{
          #----------------------------------------------------------------
          w_q <- where_q[[q]]
          der <- lapply(as.list(as.data.frame(t(w_q))),function(l){
            second_index <- matE[l[1],ifelse(l[2] == 1,2,1)]
            qprime <- second_index

            matlist_theta <- list_theta[[l[1]]]
            don <- list_Mat[[l[1]]]
            maskNA <- list_MaskNA[[l[1]]]
            if (v_distrib[l[1]] == 'bernoulli' ) { Unmdon <- (1 - don) * maskNA }
            #if (v_distrib[l[1]] %in% c('laplace','poisson','gaussian','ZIgaussian')) {Unit <- matrix(1,nrow(don),ncol(don))}
            if (v_distrib[l[1]] %in% c('laplace','poisson','gaussian','ZIgaussian')) {Unit <- maskNA}

            #if (v_distrib[l[1]] == 'ZIgaussian') { Zerosdon <- (don == 0) }
            if (v_distrib[l[1]] == 'ZIgaussian') { NonZerosdon <- (don != 0); Zerosdon <- (don == 0) * maskNA }

            #------------
            if (qprime < 1)   #if sbm
            {
              #if (v_distrib[l[1]] == 'bernoulli' ) {diag(Unmdon) <- 0 }
              #if (v_distrib[l[1]] %in% c('laplace','poisson','gaussian','ZIgaussian') ) {diag(Unit) <- 0 }
              qprime <- q
            }else{ # if lbm
              if (l[2] == 2) #functional group q at stake in rows
              {
                don <- t(don)
                if (v_distrib[l[1]] == 'bernoulli' ) {Unmdon <- t(Unmdon) }
                if (v_distrib[l[1]]  %in% c('laplace','poisson','gaussian','ZIgaussian') ) {Unit <- t(Unit) }
                #if (v_distrib[l[1]]  == 'ZIgaussian' ) {Zerosdon <- t(Zerosdon) }
                if (v_distrib[l[1]] == 'ZIgaussian') { NonZerosdon <- t(NonZerosdon); Zerosdon <- t(Zerosdon) }

                if (v_distrib[l[1]] %in% c('gaussian','ZIgaussian' )) {matlist_theta <- lapply(matlist_theta,function(theta){t(theta)})}else{matlist_theta = t(matlist_theta)}
              }
            #------------
            }

            #poisson or bernoulli likelihood
            switch(v_distrib[l[1]],
              bernoulli = {lik  = don %*% tcrossprod(tau[[qprime]],log(matlist_theta)) + Unmdon %*% tcrossprod(tau[[qprime]],log(1 - matlist_theta))},
              poisson = {lik  = -Unit %*% tcrossprod(tau[[qprime]], matlist_theta) + don %*% tcrossprod(tau[[qprime]],log(matlist_theta))},
              laplace  = {lik  = -Unit %*% tcrossprod(tau[[qprime]], log(2 * matlist_theta)) - abs(don) %*% tcrossprod(tau[[qprime]], 1/matlist_theta)},
              gaussian  = {lik = -0.5 * Unit %*% tcrossprod(tau[[qprime]], log(2 * pi * matlist_theta$var) + matlist_theta$mean^2/matlist_theta$var) - 0.5 * don^2 %*% tcrossprod(tau[[qprime]], 1/matlist_theta$var) + don %*% tcrossprod(tau[[qprime]], matlist_theta$mean/matlist_theta$var)},
              ZIgaussian  = {
                likGauss = -0.5 * NonZerosdon %*% tcrossprod(tau[[qprime]], log(2 * pi * matlist_theta$var) + matlist_theta$mean^2/matlist_theta$var) - 0.5 * ( NonZerosdon  * don^2) %*% tcrossprod(tau[[qprime]], 1/matlist_theta$var) + (NonZerosdon * don) %*% tcrossprod(tau[[qprime]], matlist_theta$mean/matlist_theta$var)
                likZeros = Zerosdon %*% tcrossprod(tau[[qprime]],log(matlist_theta$p0)) + (NonZerosdon) %*% tcrossprod(tau[[qprime]],log(1 - matlist_theta$p0))
                lik  = likGauss + likZeros}
            )



            if (second_index < 0)#sbm but non sym
            {
              don <- t(don)
              if (v_distrib[l[1]] == 'bernoulli' ) { Unmdon <- t(Unmdon) }
              if (v_distrib[l[1]]  %in% c('laplace','poisson','gaussian','ZIgaussian') ) {Unit <- t(Unit) }
              if (v_distrib[l[1]]  == 'ZIgaussian' ) {Zerosdon <- t(Zerosdon); NonZerosdon <- t(NonZerosdon) }
              if (v_distrib[l[1]] %in% c('gaussian', 'ZIgaussian' )) {
                matlist_theta <- lapply(matlist_theta,function(theta){t(theta)})
              }else{matlist_theta <- t(matlist_theta)}

              switch(v_distrib[l[1]],
                bernoulli = {lik = lik + don %*% tcrossprod(tau[[qprime]],log(matlist_theta)) + Unmdon %*% tcrossprod(tau[[qprime]],log(1 - matlist_theta))},
                poisson =   {lik = lik - Unit %*% tcrossprod(tau[[qprime]], matlist_theta)  +   don %*% tcrossprod(tau[[qprime]],log(matlist_theta))},
                laplace  =  {lik = lik - Unit %*% tcrossprod(tau[[qprime]], log(2 * matlist_theta)) - abs(don) %*% tcrossprod(tau[[qprime]], 1/matlist_theta)},
                gaussian  = {lik = lik - 0.5 * Unit %*% tcrossprod(tau[[qprime]], log(2 * pi * matlist_theta$var) + matlist_theta$mean^2/matlist_theta$var) - 0.5 * don^2 %*% tcrossprod(tau[[qprime]], 1/matlist_theta$var) + don %*% tcrossprod(tau[[qprime]], matlist_theta$mean/matlist_theta$var)},
                ZIgaussian  = {
                  likGauss = -0.5 *  NonZerosdon %*% tcrossprod(tau[[qprime]], log(2 * pi * matlist_theta$var) + matlist_theta$mean^2/matlist_theta$var) - 0.5 * ( NonZerosdon  * don^2) %*% tcrossprod(tau[[qprime]], 1/matlist_theta$var) + (NonZerosdon * don) %*% tcrossprod(tau[[qprime]], matlist_theta$mean/matlist_theta$var)
                  likZeros = Zerosdon %*% tcrossprod(tau[[qprime]],log(matlist_theta$p0)) + (NonZerosdon) %*% tcrossprod(tau[[qprime]],log(1 - matlist_theta$p0))
                #
                #
                # ZIgaussian  = {
                #   likGauss = -0.5 * ((1 - Zerosdon) * Unit) %*% tcrossprod(tau[[qprime]], log(2 * pi * matlist_theta$var) + matlist_theta$mean^2/matlist_theta$var) - 0.5 * ((1 - Zerosdon)* don^2) %*% tcrossprod(tau[[qprime]], 1/matlist_theta$var) + ((1 - Zerosdon)*don) %*% tcrossprod(tau[[qprime]], matlist_theta$mean/matlist_theta$var)
                #   likZeros = Zerosdon %*% tcrossprod(tau[[qprime]],log(matlist_theta$p0)) + (1 - Zerosdon) %*% tcrossprod(tau[[qprime]],log(1 - matlist_theta$p0))
                  lik = lik + likGauss + likZeros
                 }
              )
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


          tau[[q]] <- temp2


        }

      }
      iterVE <- iterVE + 1
      if (distTau(tau,tau_old) < valStopCrit)   stopVE <- 1

      if (iterVE == maxiterVE) {no.convergence = no.convergence + 1}; # warning(paste("Maximum number of VE iterations reached for model with v_K=", v_K,sep = ' ' ))}

    }#-------------------- END of VE Step


    #computing lik
    pseudolik <- compLikICLInt(tau,list_theta,list_pi,matE,list_Mat,n_q,v_K,v_distrib)
    vJ[iterVEM] <- pseudolik$condLik + pseudolik$margLik + pseudolik$entr
  } # ------------ end of EM var

  # convergence
  if (iterVEM == maxiterVEM) {no.convergence = no.convergence + 1}


  #computing ICL
  likICL <- compLikICLInt(tau,list_theta,list_pi,matE,list_Mat,n_q,v_K,v_distrib)

 ICL <-  likICL$condLik + likICL$margLik - 1/2 * likICL$pen

  paramEstim   <- MBMfit$new(v_K = v_K, v_distrib = v_distrib, list_pi = list_pi,list_theta = list_theta);
  names(paramEstim$list_pi) = namesFG

  names_mat <- sapply(1:nrow(dataR6$E),function(e){paste(namesFG[dataR6$E[e,1]],namesFG[dataR6$E[e,2]],sep='')})
  names(paramEstim$list_theta) = names_mat

  paramEstim$tau <- tau
  names(paramEstim$tau) = namesFG

  vJ <- vJ[1:iterVEM]
  #if (iterVEM == maxiter) { warning(paste("Maximum number of VEM iterations reached for model with v_K=", v_K,sep = ' ' ))}
  return( list(paramEstim = paramEstim,ICL = ICL,vJ = vJ, convergence = (no.convergence == 0)))
}

