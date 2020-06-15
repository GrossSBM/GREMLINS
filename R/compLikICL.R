#' compute the Integrated likeilhood and the ICL criteria for the MBM
#'
#' @param paramEstim Estimated parameters of MBM
#' @param list_Net A list of network
#' @param v_distrib  Type of proababilistic distributions in each network : if 0/1 then Bernoulli, if counting then Poisson. My default  = Bernoulli.
#'                   Must give a vector whose length is the number of networks in list_Net
#' @return Pseudo-Likelihood and penalty
#' @examples
#'
#' @export


compLikICL  <- function(paramEstim,list_Net,v_distrib = NULL)
{

  dataR6 = formattingData(list_Net,v_distrib)

  list_MaskNA <- lapply(dataR6$mats,function(m){1 * (1 - is.na(m))})

  tau <- paramEstim$tau
  list_pi <- paramEstim$list_pi
  if (is.null(list_pi)) { list_pi <- lapply(tau,colMeans) }

  list_theta <-  paramEstim$list_theta

  if (is.null(v_distrib)) {v_distrib <- paramEstim$v_distrib}

  res <- compLikICLInt(tau,list_theta,list_pi,dataR6$Ecode,dataR6$mats,list_MaskNA,dataR6$v_NQ,paramEstim$v_K,v_distrib)



  return(res)
}
