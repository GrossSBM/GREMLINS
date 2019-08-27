#' Compare two classifications on all the  Functional groups
#'
#' @param  classif1 : list a length n_FG.
#' @param  classif2 : list a length n_FG.
#' @return Adjusted Rand Index (ARI) for each Functional Group.
#' @examples
#' nFG <- 3;
#' vK <- c(4,5,2) ;
#' vNQ <- c(100,40,50);
#' classif1 <- lapply(1:nFG,function(q){sample(1:vK[q],vNQ[q],replace=TRUE)})
#' classif2 <- classif1
#' classif2[[2]] <-  sample(1:vK[2],vNQ[2],replace=TRUE)
#' resCompar <- comparClassif (classif1,classif2)
#' @export
#'
comparClassif <- function(classif1,classif2){

  if (length(classif1) != length(classif2) )
  {
    stop('classification can not be compared because objects of different sizes')
  }
  ARI <- vapply(1:length(classif1),function(q){adjustedRandIndex(classif1[[q]],classif2[[q]])},1)
  return(ARI)
}
