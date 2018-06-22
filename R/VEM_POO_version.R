
#initialization provides a first classification and corresponding tau
#method can be random, CAH, given in which case add given classif


VEM_gen_BM <- function(dataR6,classif.init,tau.init=NULL)
  #data :  coll_interaction type object
  #classif  : liste de classifications (Z) au sein des groupes fonctionnels fg (de longueur Q )

{
  #checking the dimensions of matrices and extracting number of indiviuals
  #for a given q functional groups, gives the list of matrix in row or in columns where it plays a role

  #browser()
  where_q <- dataR6$where
  n_q <- dataR6$v_NQ
  cardE <- dataR6$card_E
  mat_E <- dataR6$Ecode ### a remettre en public?????
  list_Mat <- dataR6$mats
  vdistrib <- dataR6$vdistrib

  #some constants
  eps <- 2*.Machine$double.eps
  #eps <- 0.001
  val_stopcrit <- 1e-6
  val_stopcrit <- 1e-10



  ##  initialisation
  vK <- calc_vK(classif.init)
  tau <- tau.init
  if (is.null(tau)) {
    tau <-  lapply(1:dataR6$Q,function(j){
      mat <- matrix(eps,n_q[j],vK[j],byrow = TRUE)
      mat[cbind(1:n_q[j],classif.init[[j]])] <- 1 - eps
      #normalize tau
      mat <- mat/rowSums(mat)
      return(mat)
    })
  }




  #for computations when SBM putting 0 on the diagonal
  indSBM <- which(mat_E[,2] < 1)
  for (i in indSBM) {diag(list_Mat[[i]])  <- 0}


  #entering VEM
  maxiter <- 1000
  stopcrit <- 0
  iter_VEM <- 0

  vJ  <- numeric(maxiter)

  #for stopping criterion
  ltheta <- lapply(1:cardE,function(s)
  {
    gr <- mat_E[s,1]
    gc <- mat_E[s,2]
    if (gc < 1) gc <- gr
    return(matrix(Inf,vK[gr],vK[gc]))
  })


  ######################
  # Algo begins
  #####################
  while (iter_VEM < maxiter & stopcrit == 0)
  {
    iter_VEM <- iter_VEM + 1
    #if(iter_VEM%%100==0){print(paste("Iteration of VEM",iter_VEM,sep=' : '))}


    #--------------------------------   M step
    ltheta_old <- ltheta

    lpi = lapply(tau,colMeans)
    ltheta  = lapply(1:cardE,function(j){
      gr <- mat_E[j,1]
      gc <- mat_E[j,2]


      if (gc < 1) {  #for sbm sym or notsym
        gc <- gr
        #useful matrix
        Unitmdiag <- matrix(1,nrow = n_q[gr],ncol = n_q[gc])
        diag(Unitmdiag) <- 0
        #bernoulli or poisson distribution same expression for M step
        lthetac = t(tau[[gr]]) %*% list_Mat[[j]] %*% tau[[gc]] / (t(tau[[gr]]) %*% (Unitmdiag) %*% tau[[gc]])
      }
      else #for lbm
      {
        Unit <- matrix(1,nrow = n_q[gr],ncol = n_q[gc])
        lthetac <-  t(tau[[gr]]) %*% list_Mat[[j]] %*% tau[[gc]] / (t(tau[[gr]]) %*% (Unit) %*% tau[[gc]])
      }
      return(lthetac)})

    #prevent the values from being too close from 0 or 1
    lpi <- lapply(lpi,readjust_pi,eps)
    ltheta <- lapply(ltheta,readjust_theta,eps)


    #stop criterion
    if (distltheta(ltheta,ltheta_old) < val_stopcrit) stopcrit <- 1


    #--------------------------------   VE step : boucle
    iterVE=0
    stopVE=0
    while (iterVE < maxiter&stopVE == 0)
    {
      #VE step
      tau_old=tau #useful ?
      for (q in 1:dataR6$Q)
      {
        w_q=where_q[[q]]
        der=lapply(as.list(as.data.frame(t(w_q))),function(l)
        {
          second_index=mat_E[l[1],ifelse(l[2]==1,2,1)]
          qprime=second_index

          don=list_Mat[[l[1]]] #pb a regler ligne ou colonne et si sbm enlever de 1-don la diag
          Unmdon=1-don
          matltheta=ltheta[[l[1]]]

          if (qprime<1)   #if sbm
          {
            diag(Unmdon)=0
            qprime=q
          }
          else  # if lbm
          {
            if (l[2]==2) #functional group q at stake in rows
            {
              don=t(don)
              Unmdon=t(Unmdon)
              matltheta=t(matltheta)
            }
          }

          #poisson or bernoulli likelihood
          switch(vdistrib[l[1]],
            bernoulli={lik=don%*%tau[[qprime]]%*%t(log(matltheta))+Unmdon%*%tau[[qprime]]%*%t(log(1-matltheta))},
            poisson={stop("Codes non ecrits pour les lois de poisson")}
          )



          if (second_index < 0)#sbm but non sym
          {
            don <- t(don)
            Unmdon <- t(Unmdon)
            matltheta <- t(matltheta)
            switch(vdistrib[l[1]],
              bernoulli = {lik = lik + don %*% tau[[qprime]] %*% t(log(matltheta)) + Unmdon %*% tau[[qprime]] %*% t(log(1 - matltheta))},
              poisson = {stop("Codes non ecrits pour les lois de poisson")})
          }
          return(lik)
        })
        L <- (Reduce('+',der))

        B <- L + matrix(log(lpi[[q]]),nrow = nrow(tau[[q]]),ncol = vK[q],byrow = TRUE)
        B <- B - max(B)

        temp <- exp(B)
        temp2 <- temp/rowSums(temp)
        temp2[temp2 < eps] <- eps
        temp2[temp2 > 1 - eps] <- 1 - eps
        temp2 = temp2/rowSums(temp2)
        tau[[q]] = temp2
      }

      #boucle VE
      iterVE = iterVE + 1
      if (disttau(tau,tau_old) < val_stopcrit)   stopVE <- 1
    }
    #

    #computing lik
    pseudolik <- comp_lik_ICL(tau,ltheta,lpi,mat_E,list_Mat,n_q,vK)
    vJ[iter_VEM] <- pseudolik$condlik + pseudolik$marglik + pseudolik$entr
  }

  #computing ICL
  likicl <- comp_lik_ICL(tau,ltheta,lpi,mat_E,list_Mat,n_q,vK)

  icl <-  likicl$condlik + likicl$marglik - 1/2*likicl$pen

  param_estim   <- genBMfit$new(vK = vK, vdistrib = vdistrib, lpi = lpi,ltheta = ltheta);
  param_estim$tau <- tau
  vJ <- vJ[1:iter_VEM]
  return(list(param_estim = param_estim,ICL = icl,vJ = vJ))
}

