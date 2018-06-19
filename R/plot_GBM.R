#' Plot the summary of the estimated generalized blockmodel
#'
#' @param listNet A list of network (defined via the function DefineNetwork)
#' @param fitted_GBM A fitted Generalized BlockModel
#' @param mycol A list of colors for the functional groups
#' @param thres A threshold under which edges correponding to probability of connections are not plotted
#' @examples
#' npc1 <- 20 # nodes per class
#' Q1 <- 3 # classes
#' n1 <- npc1 * Q1 # nodes
#' Z1 <- diag(Q1)%x%matrix(1,npc1,1)
#' P1 <- matrix(runif(Q1*Q1),Q1,Q1)
#' A <- 1*(matrix(runif(n1*n1),n1,n1)<Z1%*%P1%*%t(Z1)) ## adjacency matrix
#' Agr <- DefineNetwork(A,"diradj","FG1","FG1")
#' npc2 <- 30 # nodes per class
#' Q2 <- 2 # classes
#' n2 <- npc2 * Q2 # nodes
#' Z2 <- diag(Q2)%x%matrix(1,npc2,1)
#' P2 <- matrix(runif(Q1*Q2),Q1,Q2)
#' B <- 1*(matrix(runif(n1*n2),n1,n2)<Z1%*%P2%*%t(Z2)) ## adjacency matrix
#' Bgr <- DefineNetwork(B,"inc","FG1","FG2")
#' res <- MultipartiteBM(list(Agr,Bgr),namesFG = NULL,vKmin = 1,vKmax = 10,vKinit = NULL,verbose = TRUE, save=FALSE)
#' plot_GBM(list(Agr,Bgr),res$param_estim)
#' @export

plot_GBM = function(listNet,fitted_GBM,mycol=NULL,thres=0.01){

  dataR6 = FormattingData(listNet)
  Q <-  dataR6$Q
  vK_estim <- fitted_GBM$vK
  type_node <- lapply(1:Q,function(q){rep(dataR6$namesfg[q],vK_estim[q])})
  label_node <- lapply(1:Q,function(q){1:vK_estim[q]})
  if(is.null(mycol)){mycol <-  palette("default"); mycol <- mycol[-1]}
  col_node <- lapply(1:Q,function(q){rep(mycol[q],vK_estim[q])})
  size_node <- lapply(1:Q,function(q){fitted_GBM$lpi[[q]]})
  cum_vK =c(0,cumsum(vK_estim))
  code_node <- lapply(2:(Q+1),function(q){seq(cum_vK[q-1]+1,cum_vK[q],1)})


  N = sum(vK_estim)
  list_edges <- lapply(1:nrow(dataR6$E),function(i){
    q.row = dataR6$E[i,1]
    q.col = dataR6$E[i,2]
    ltheta_i <- fitted_GBM$ltheta[[i]]
    c1 <- rep(code_node[[q.row]],times=vK_estim[q.col])
    c2 <- rep(code_node[[q.col]],each = vK_estim[q.row])
    edges_i <- cbind(c1,c2,c(ltheta_i))
    edges_i <- as.data.frame(edges_i)
    edges_i$type <- rep(dataR6$type_inter[i],length(c1))
    return(edges_i)})
  all_edges <- do.call("rbind", list_edges)

  # all_edges <- rbind(all_edges,c(2,2,10))
  # all_edges$type[nrow(all_edges)] = "diradj"
  w <- which(all_edges[,3]>thres)
  edges <- all_edges[w,c(1,2)]
  curved <- 0.3*(all_edges[w,4]=="diradj")
  curved <- runif(length(all_edges[w,4] == "diradj"),0,3)*(all_edges[w,4]=="diradj")




  if(sum(dataR6$type_inter=="diradj")==0){
    G<- make_empty_graph() + vertices(unlist(code_node))
    G <- G %>%set_vertex_attr("label",value=unlist(label_node)) %>%set_vertex_attr("color",value=unlist(col_node))
    G <- G%>%set_vertex_attr("size",value=sqrt(unlist(size_node))*20+4)
    G <- G%>%set_edge_attr("width",value=all_edges[all_edges[,3]>=thres,3] )
    G <- G + edges(c(t(edges)))
    G <- as.undirected(G)

  }
  if(sum(dataR6$type_inter=="diradj")>0){
    G<- make_empty_graph() + vertices(unlist(code_node))
    G <- G %>%set_vertex_attr("label",value=unlist(label_node)) %>%set_vertex_attr("color",value=unlist(col_node))
    G <- G%>%set_vertex_attr("size",value=sqrt(unlist(size_node))*20+4)
    G <- G%>%set_edge_attr("width",value=all_edges[all_edges[,3]>=thres,3] )
    # G <- G%>%set_edge_attr("curved", E(.)[1%--%3,3%--%1],value = 0.4)
    G <- G + edges(c(t(edges)))
  }
  m <- max(all_edges[,3])

  plot(G,edge.width=all_edges[all_edges[,3]>thres,3]/m*5,edge.curved=curved)




  #coords <- layout_in_circle(G)

  legend("left", c(dataR6$namesfg), col=mycol[1:Q], border = "black", lty=1, lwd=4)
}

