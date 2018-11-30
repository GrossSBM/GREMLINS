#' compute the Integrated likeilhood and the ICL criteria for the MBM
#'
#' @param param_est Estimated parameters of MBM
#' @param list_Net A list of network
#' @return Pseudo-Likelihood, penalty
#' @examples
#' A <- matrix(rbinom(100,1,.2),10,10)
#' type <- "diradj"
#' defineNetwork(A,"diradj","FG1","FG1")
#' @export


compLikICLExport  <- function(paramEstim,list_Net)
{

  dataR6 = formattingData(list_Net)
  matE <- dataR6$Ecode
  n_q <- dataR6$v_NQ
  cardE <- dataR6$cardE
  list_Mat <- dataR6$mats

  #browser()
  v_K <- paramEstim$v_K
  tau <- paramEstim$tau
  piEstim <- paramEstim$list_pi
  thetaEstim <- paramEstim$list_theta

  if (is.null(piEstim)) { piEstim <- lapply(tau,colMeans) }
  if (is.null(thetaEstim)) {
    thetaEstim  <- lapply(1:cardE,function(j){
      gr <- matE[j,1]
      gc <- matE[j,2]


      if (gc < 1) {  #for sbm sym or notsym
        gc <- gr
        #useful matrix
        Unitmdiag <- matrix(1,nrow = n_q[gr],ncol = n_q[gc])
        diag(Unitmdiag) <- 0
        #bernoulli or poisson distribution same expression for M step
        list_thetac <- t(tau[[gr]]) %*% list_Mat[[j]] %*% tau[[gc]] / (t(tau[[gr]]) %*% (Unitmdiag) %*% tau[[gc]])
      }
      else #for lbm
      {
        Unit <- matrix(1,nrow = n_q[gr],ncol = n_q[gc])
        list_thetac <-  t(tau[[gr]]) %*% list_Mat[[j]] %*% tau[[gc]] / (t(tau[[gr]]) %*% (Unit) %*% tau[[gc]])
      }
      return(list_thetac)})
  }

  #  E_Y[log l(X | Z)]
  condLiks = sapply(1:cardE,function(e)
  {
    gr <- matE[e,1]  # FG on row
    gc <- matE[e,2]  # FG on column
    don <- list_Mat[[e]]
    Unmdon <- 1 - don

    #-----
    facteur = 1
    if (gc < 1) # then sbm
    {
      if (gc == 0)  facteur = 1/2 #sbm sym # else -1 sbm non sym
      gc = gr
      diag(Unmdon) = 0
    }
    #-----
    prov = (tau[[gr]]) %*% log(thetaEstim[[e]]) %*% t(tau[[gc]])
    prov1m = (tau[[gr]]) %*% log(1 - thetaEstim[[e]]) %*% t(tau[[gc]])
    return((sum(don * prov) + sum(Unmdon*prov1m)) * facteur)
  }
  )
  condLik = sum(condLiks)

  #-----  E_Y[log l(Z)]
  margLiks = sapply(1:dataR6$Q,function(q)
  {
    return(sum(tau[[q]]*matrix(log(piEstim[[q]]),nrow(tau[[q]]),ncol(tau[[q]]),byrow = TRUE)))
  })
  margLik = sum(margLiks)

  #----- Entropy

  entrs = sapply(1:dataR6$Q,function(q)
  {
    return(-sum(log(tau[[q]])*tau[[q]]))
  })
  entr = sum(entrs)

  #------ Penalty term
  penMats = vapply(1:cardE,function(s)
  {
    gr = matE[s,1]
    gc = matE[s,2]
    if (gc > 0) #LBM
    {
      #return(c(v_K[gr]*v_K[gc], prod(dim(list_Mat[[s]]))))
      return(c(v_K[gr]*v_K[gc], n_q[gr] * n_q[gc]))
    }
    if (gc == 0) #SBM sym
    {
      #return(c(v_K[gr]*(v_K[gr]+1)/2, nrow(list_Mat[[s]])*(nrow(list_Mat[[s]])-1)/2))
      #
      return(c(v_K[gr]*(v_K[gr] + 1) * 0.5, 0.5 * n_q[gr]*(n_q[gr] - 1)))

    }
    if (gc == -1) #SBM non sym
    {
      return(c(v_K[gr] * (v_K[gr]), n_q[gr] * (n_q[gr] - 1)))
    }

  },c(1.0,2.0)
  )


  penICL = sum((v_K - 1) * log(n_q))  + sum(penMats[1,]) * log(sum(penMats[2,]))

  return(list(condLik = condLik,margLik = margLik,entr = entr,pen = penICL))
}
