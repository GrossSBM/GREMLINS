#' Define a network providing its matrix of interactions and specifying the functions groups in row and col.
#'
#' @param mat An adjacency matrix (symmetric or not) or an incidence matrix
#' @param typeInter Type of the matrix, choice between "inc" (incidence), "adj" (adjacency) and "diradj" (directed adjacency)
#' @param rowFG Name of the functional group in row
#' @param colFG Name of the function group in column
#' @return a list object formatted for the GREMLIN package
#' @examples
#' A <- matrix(rbinom(100,1,.2),10,10)
#' type <- "diradj"
#' defineNetwork(A,"diradj","FG1","FG1")
#' @export

defineNetwork = function(mat,typeInter,rowFG,colFG)
{
  if (typeInter %in% c("inc","adj","diradj") == F) {stop("not allowed type")}
  if (typeInter == "adj" & !isSymmetric(mat)) {stop("not symmetric adjacency matrix")}
  else return(list(mat = mat,typeInter = typeInter,rowFG = rowFG,colFG = colFG))
}

