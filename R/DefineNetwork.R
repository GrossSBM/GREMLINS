#' Define a network
#'
#' @param mat An adjacency matrix (symmetric or not) or an incidence matrix
#' @param type Type of the matrix, choice between "inc" (incidence), "adj" (adjacency) and "diradj" (directed adjacency)
#' @param rowFG Name of the functional group in row
#' @param colFG Name of the function group in column
#' @return a list object formatted for the GREMLIN package
#' @examples
#' A <- matrix(rbinom(100,1,.2),10,10)
#' type <- "diradj"
#' DefineNetwork(A,"diradj","FG1","FG1")
#' @export

DefineNetwork = function(mat,type,rowFG,colFG)
{
  if (type %in% c("inc","adj","diradj") == F) {stop("not allowed type")}
  if (type == "adj" & !isSymmetric(mat)) {stop("not symmetric adjacency matrix")}
  else return(list(mat = mat,type = type,rowFG = rowFG,colFG = colFG))
}

