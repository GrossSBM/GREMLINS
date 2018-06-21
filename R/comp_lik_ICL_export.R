#' compute the Integrated likeilhood and the ICL criteria for the MBM
#'
#' @param tau Posterior probabilities of clustering
#' @param pi Type of the matrix, choice between "inc", "adj" and "diradj"
#' @param alpha Name of the functional group in row
#' @param colFG Name of the function group in column
#' @return a list object formatted for the GREMLIN package
#' @examples
#' A <- matrix(rbinom(100,1,.2),10,10)
#' type <- "diradj"
#' DefineNetwork(A,"diradj","FG1","FG1")
#' @export




#computing ICL and likelihood
comp_lik_ICL_export  <- function(res_estim,listNet)
{



  dataR6 = FormattingData(listNet)

  #dataR6 <-  FormattingData(listNet)

  mat_E <- dataR6$Ecode
  n_q <- dataR6$v_NQ
  cardE <- dataR6$card_E

  list_Mat <- dataR6$mats


  vK <- res_estim$param_estim$vK
  tau <- res_estim$param_estim$tau
  pi_estim <- res_estim$param_estim$lpi
  theta_estim <- res_estim$param_estim$ltheta

  #  E_Y[log l(X | Z)]
  condliks = sapply(1:cardE,function(e)
  {
    gr <- mat_E[e,1]  # FG on row
    gc <- mat_E[e,2]  # FG on column
    don <- list_Mat[[e]]
    Unmdon <- 1 - don

    #-----
    facteur = 1
    if (gc<1) # then sbm
    {
      if (gc == 0)  facteur=1/2 #sbm sym # else -1 sbm non sym
      gc = gr
      diag(Unmdon) = 0
    }
    #-----
    prov = (tau[[gr]]) %*% log(theta_estim[[e]]) %*% t(tau[[gc]])
    prov1m = (tau[[gr]])%*%log(1-theta_estim[[e]])%*%t(tau[[gc]])
    return((sum(don*prov) + sum(Unmdon*prov1m))*facteur)
  }
  )
  condlik=sum(condliks)

  #-----  E_Y[log l(Z)]
  likmargs=sapply(1:dataR6$Q,function(q)
  {
    return(sum(tau[[q]]*matrix(log(pi_estim[[q]]),nrow(tau[[q]]),ncol(tau[[q]]),byrow = TRUE)))
  })
  likmarg=sum(likmargs)

  #----- Entropy

  entros=sapply(1:dataR6$Q,function(q)
  {
    return(-sum(log(tau[[q]])*tau[[q]]))
  })
  entro = sum(entros)

  #------ Penalty term
  pen_mats=vapply(1:cardE,function(s)
  {
    gr=mat_E[s,1]
    gc=mat_E[s,2]
    if (gc>0) #LBM
    {
      #return(c(vK[gr]*vK[gc], prod(dim(list_Mat[[s]]))))
      return(c(vK[gr]*vK[gc], n_q[gr] * n_q[gc]))
    }
    if (gc==0) #SBM sym
    {
      #return(c(vK[gr]*(vK[gr]+1)/2, nrow(list_Mat[[s]])*(nrow(list_Mat[[s]])-1)/2))
      #
      return(c(vK[gr]*(vK[gr] + 1) * 0.5, 0.5 * n_q[gr]*(n_q[gr] - 1)))

    }
    if (gc==-1) #SBM non sym
    {
      return(c(vK[gr]*(vK[gr]), n_q[gr]*(n_q[gr]-1)))
    }

  },c(1.0,2.0)
  )


  penICL=sum((vK-1)*log(n_q))  + sum(pen_mats[1,])*log(sum(pen_mats[2,]))


  return(list(condlik=condlik,marglik=likmarg,entr=entro,pen=penICL))

}
