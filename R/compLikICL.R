#' compute the Integrated likeilhood and the ICL criteria for the MBM
#'
#' @param param_est Estimated parameters of MBM
#' @param list_Net A list of network
#' @return Pseudo-Likelihood, penalty
#' @examples
#'
#' @export


compLikICL  <- function(paramEstim,list_Net,v_distrib = NULL)
{

  dataR6 = formattingData(list_Net,v_distrib)



  tau <- paramEstim$tau
  list_pi <- paramEstim$list_pi
  if (is.null(list_pi)) { list_pi <- lapply(tau,colMeans) }

  list_theta <-  paramEstim$list_theta

  if (is.null(v_distrib)) {v_distrib <- paramEstim$v_distrib}

  res <- compLikICLInt(tau,list_theta,list_pi,dataR6$Ecode,dataR6$mats,dataR6$v_NQ,paramEstim$v_K,v_distrib)



  return(res)
}
