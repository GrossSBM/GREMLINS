#' Define a network providing its matrix of interactions and specifying the functions groups in row and col.
#'
#' @param mat An adjacency matrix (symmetric or not) or an incidence matrix
#' @param typeInter Type of the matrix, choice between "inc" (incidence), "adj" (adjacency) and "diradj" (directed adjacency)
#' @param rowFG Name of the functional group in row
#' @param colFG Name of the function group in column
#' @param selfLoops Only used for adjacency matrices ("adj" or "diradj"). If \code{FALSE} (default),
#'               the diagonal of \code{mat} is set to \code{NA}: self-loops are considered as non observed and
#'               do not contribute to the inference (as in the packages \pkg{blockmodels} and \pkg{sbm}).
#'               If \code{TRUE}, the diagonal is kept and modelled as any other entry.
#' @return a list object formatted for the GREMLINS package
#' @examples
#' A <- matrix(rbinom(100,1,.2),10,10)
#' type <- "diradj"
#' defineNetwork(A,"diradj","FG1","FG1")
#' @export

defineNetwork = function(mat,typeInter,rowFG,colFG,selfLoops = FALSE)
{
  if (typeInter %in% c("inc","adj","diradj") == F) {stop("not allowed type")}
  if (typeInter == "adj" & !isSymmetric(mat)) {stop("not symmetric adjacency matrix")}
  if (typeInter != "inc" & !selfLoops) {diag(mat) <- NA}
  return(list(mat = mat,typeInter = typeInter,rowFG = rowFG,colFG = colFG,selfLoops = selfLoops))
}

