#' Plot the mesoscopic view  of the estimated MBM
#'
#' Plot the mesocopic view of the multipartite network obtained by the Genreliazed block models. resMBM is the results of of MBM fitting (output of multipartiteBMFixedModel for given numbers of clusters or multipartiteBM if the number of clusters is selected)
#'
#' @param resMBM A fitted Generalized BlockModel
#' @param whichModel The index corresponding to the model to plot (default is 1, the best model)
#' @param mycol A list of colors for the functional groups
#' @param thres A threshold under which edges correponding to probability of connections are not plotted
#' @param cex graphical parameter:
#' @param maxCurved graphical parameter : curvature of the edges
#'
#' @examples
#' v_K <- c(3,2,2)
#' n_FG <- 3
#' list_pi <- vector("list", 3);
#' list_pi[[1]] <- c(0.4,0.3,0.3); list_pi[[2]] <- c(0.6,0.4); list_pi[[3]]  <- c(0.6,0.4)
#' E  = rbind(c(1,2),c(2,3),c(2,2))
#' v_distrib <- c('bernoulli','poisson','poisson')
#' typeInter <- c( "inc", "inc"  ,  "adj" )
#' list_theta <- list()
#' list_theta[[1]] <- matrix(rbeta(v_K[E[1,1]] * v_K[E[1,2]],1.5,1.5 ),nrow = v_K[E[1,1]], ncol = v_K[E[1,2]])
#' list_theta[[2]] <- matrix(rgamma(v_K[E[2,1]] * v_K[E[2,2]],7.5,1 ),nrow = v_K[E[2,1]], ncol = v_K[E[2,2]])
#' list_theta[[3]] <- matrix(rgamma(v_K[E[3,1]] * v_K[E[3,2]],7.5,1 ),nrow = v_K[E[3,1]], ncol = v_K[E[3,2]])
#' list_theta[[3]] <- 0.5*(list_theta[[3]] + t(list_theta[[3]])) # symetrisation for network 3
#' v_NQ = c(100,50,40)
#' dataSim <-  rMBM(v_NQ ,E , typeInter, v_distrib, list_pi, list_theta, seed=NULL, namesFG= c('A','B','D'),keepClassif = FALSE)
#' list_Net <- dataSim$list_Net
#' resMBM <- multipartiteBMFixedModel(list_Net,namesFG = c('A','B','D'), v_K = c(3,2,2),v_distrib = v_distrib)
#' plotMBM(resMBM,mycol=c('magenta','cyan','yellow'))
#' @export

plotMBM = function(resMBM,whichModel = 1, mycol = NULL, thres = 0.01,  cex=1, maxCurved=3){

  list_Net <- resMBM$list_Net
  Q <- length(resMBM$fittedModel[[1]]$paramEstim$list_pi)

  dataR6 <- formattingData(list_Net,v_distrib = resMBM$fittedModel[[whichModel]]$paramEstim$v_distrib)
  labelFG <- substr(dataR6$namesFG,1,1)
  i = 1;
  while ( length(unique(labelFG))<Q & (i <= Q)){i = i + 1; labelFG <- substr(dataR6$namesFG,1,i)}

  if ((length(resMBM$fittedModel) == 1) & (whichModel > 1)) {stop('Only one  fitted model in resMBM. You can not select a whichModel not equal to 1')}



  nbNet <- length(resMBM$fittedModel[[1]]$paramEstim$list_theta)
  param <- resMBM$fittedModel[[whichModel]]$paramEstim
  v_K <- param$v_K

  labelNode <- lapply(1:Q,function(q){paste(labelFG[q],1:v_K[q],sep='.')})
  if (is.null(mycol)) {mycol <-  palette("default");  mycol <- mycol[-1]}
  colNode <- lapply(1:Q,function(q){rep(mycol[q],v_K[q])})
  sizeNode <- lapply(1:Q,function(q){param$list_pi[[q]]})
  cumVK <-  c(0,cumsum(v_K))
  codeNode <- lapply(2:(Q + 1),function(q){seq(cumVK[q - 1] + 1,cumVK[q],1)})



  list_edges <- lapply(1:nbNet,function(i){
    qRow <- dataR6$E[i,1]
    qCol <- dataR6$E[i,2]
    list_theta_i <- param$list_theta[[i]]
    c1 <- rep(codeNode[[qRow]],times = v_K[qCol])
    c2 <- rep(codeNode[[qCol]],each = v_K[qRow])
    edges_i <- cbind(c1,c2,c(list_theta_i))
    edges_i <- as.data.frame(edges_i)
    edges_i$type <- rep(dataR6$typeInter[i],length(c1))
    return(edges_i)})
  allEdges <- do.call("rbind", list_edges)


  allEdges$arrow_mode <- rep(0,nrow(allEdges))  # directed or nont directed
  allEdges$arrow_mode[allEdges$type == "diradj"] = 2

  # allEdges <- rbind(allEdges,c(2,2,10))
  # allEdges$type[nrow(allEdges)] = "diradj"
  w <- which(allEdges[,3] > thres)
  edges <- allEdges[w,c(1,2)]
  curved <- 0 * (allEdges[w,4] == "diradj")

  curved <- runif(length(allEdges[w,4] == "diradj"),0,maxCurved)*(allEdges[w,4] == "diradj")



  G <- make_empty_graph() + vertices(unlist(codeNode))
  V(G)$label.cex = 1

  G <- G  %>% set_vertex_attr("label",value = unlist(labelNode))
  G <- G  %>% set_vertex_attr("color",value = unlist(colNode))
  G <- G  %>% set_vertex_attr("size",value = sqrt(unlist(sizeNode)) * 40 + 2)
  G <- G  %>% set_edge_attr("width",value = allEdges[allEdges[,3] >= thres,3] )
  G <- G + edges(c(t(edges)))
  G <- G  %>% set_graph_attr("layout" , layout_with_lgl)
  m <- max(allEdges[,3])
  plot(G,edge.width = allEdges[allEdges[,3] > thres,3] / m * 5 ,edge.curved = curved, edge.arrow.mode = allEdges$arrow_mode[w])
  legend("topleft", c(dataR6$namesFG), col = mycol[1:Q], border = "black", lty = 1, lwd = 4,cex = cex)


}

