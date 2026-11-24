#' Adjusting an extended SBM to Multipartite networks
#'
#' Generalized multipartite networks  consist in the joint observation of several networks implying some common pre-specified groups of individuals. GREMLINS adjusts an adapted version of  the popular stochastic block model to multipartite networks, as described in Bar-Hen, Barbillon and Donnet (2022).
#'
#' **GREMLINS is the inference engine of the package `sbm`. We recommend using it through `sbm::estimateMultipartiteSBM()`**,
#' which implements the same model and algorithm with a simpler interface and plot functions
#' (see `vignette("Introduction", package = "GREMLINS")` for the correspondence between the functions).
#'
#' The GREMLINS package provides the following top-level  major functions:
#'\describe{
#'   \item{\code{\link{defineNetwork}}}{ a function to define carefully a single network.}
#'   \item{\code{\link{rMBM}}}{ a function to simulate a collection of networks involving common functional groups of entities (with various emission distributions).}
#'   \item{\code{\link{multipartiteBM}}}{ a function to perform inference (model selection and estimation ) of SBM for a multipartite network.}
#'   \item{\code{\link{multipartiteBMFixedModel}}}{ a function to estimate the parameters of SBM for a multipartite network for fixed numbers of blocks}
#'}
#'
#' We also provide some additional functions useful to analyze the results:
#' \describe{
#'   \item{\code{\link{extractClustersMBM}}}{ a function to extract the clusters in each functional group}
#'   \item{\code{\link{comparClassif}}}{ a function to compute the Adjusted Rand Index (ARI) between two classifications}
#'   \item{\code{\link{predictMBM}}}{ a function to compute the predictions once the model has been fitted}
#'   \item{\code{\link{compLikICL}}}{ a function to	compute the Integrated Likelihood and the ICL criteria for the MBM}
#' }
#'
#' @name GREMLINS
#' @keywords internal
"_PACKAGE"
#' @author Pierre Barbillon, Sophie Donnet
#' @references Bar-Hen, A., Barbillon, P. and Donnet, S. (2022). Block models for generalized multipartite networks: Applications in ecology and ethnobiology. Statistical Modelling, 22(4), 273-296. \doi{10.1177/1471082X20963254}
#' @importFrom R6 R6Class
#' @import progress igraph blockmodels aricode
#' @name GREMLINS
NULL

