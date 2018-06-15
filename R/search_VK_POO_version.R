search_KQ <- function(data,vKinit,Kmin=NULL,Kmax=NULL,nb_cores=NULL,verbose = TRUE){
  #   
  
  
   
  os <- Sys.info()["sysname"]
  if ((os != 'Windows') & (is.null(nb_cores))) {nb_cores <- detectCores(all.tests = FALSE, logical = TRUE) %/% 2}
  vdistrib <- data$vdistrib
  
  #------------------------  Bouds of the numer of clusters 
  if (is.null(Kmin)) {Kmin = 1; print("The minimum number of clusters has been set to 1")}
  if (is.null(Kmax)) {Kmax = 10; print("The maximum number of clusters has been set to 10")}
  
  if (length(Kmin) == 1) {Kmin = rep(Kmin,data$Q)} else {if (length(Kmin) != data$Q) {stop("Lower bounds on vK are not of the adequate size")}}
  if (length(Kmax) == 1) {Kmax = rep(Kmax,data$Q)} else {if (length(Kmax) != data$Q) {stop("Upper bounds on vK are not of the adequate size")}}
  
  
  #------ 

  if (length(vKinit) == 1)
  {
    vKinit_temp <- switch(vKinit,
           min =  Kmin,
           mean = floor((Kmax + Kmin)/2),
           vKinit)
     vKinit = vKinit_temp
  }
  
  
  # 
  
  if (length(vKinit) != data$Q) { stop('Length of vKinit incorrect') }
  
  if (verbose) { print(paste("Searching the numbers of blocs starting from",paste(as.character(vKinit),collapse = " "),"clusters",sep = " "))}
  #------------------------  Initialisation of the number of clusters and the classification. 
  
  param.init <- genBMfit$new(vK = vKinit,vdistrib = vdistrib)

  #if (is.null(vKinit)) {if(length(Kmin == 1)){param.init$vK = rep(Kmin,data$Q)} else {param.init$vK = Kmin}} else{param.init$vK = vKinit}
  classif.init = initialize(data,param.init,method="CAH")$groups
  
  #----------------------   Initialisation of the algorithm
  
  niter.search <- 0; 
  ICL.c <- -Inf;
  classif.c <- classif.init;
  classif.new <- classif.c
  estim.new <- data$estime(classif.init);
  param.new <- estim.new$param_estim
  classif.new <- lapply(1:data$Q,function(q){max.col(param.new$tau[[q]])})
  ICL.new <- estim.new$ICL
  
  if (verbose) {
    mess <- paste(round(c(calc_vK(classif.new))),collapse = " " )
    print(paste("ICL :",round(ICL.new,2),". Nb of clusters: ", mess,sep = " "))
  }
  
  vec.ICL = ICL.new
  RES = list()
  RES[[niter.search + 1]] <- estim.new;  
  #----------------------------------------------  ALGORITHM
  while((ICL.new > ICL.c)&(niter.search<1000)){
    #
    niter.search <- niter.search + 1;
    
    #
    ICL.c <- ICL.new
    classif.c <- classif.new
    
    # list of new classif deriving from the plitting of one cluster or the merging of 2 clusters in cluterisation classif.c 
    # (These clusterisations will serve as initisalition of the EM algorithm for models) 
    list_classif_init <- sequential_initialize(classif.c ,data,Kmin,Kmax,os); 
    L = length(list_classif_init)
    
    
    if (os == "Windows") {
      all_estim <- lapply(1:L,function(l){
        print(l);
        estim.c.l <- data$estime(list_classif_init[[l]])})
    }else{
      all_estim <- mclapply(1:L,function(l){estim.c.l <- data$estime(list_classif_init[[l]])},mc.cores = nb_cores)
    }
    
    ICL.vec <- vapply(1:L,function(l){all_estim[[l]]$ICL},1)
    ICL.new <- max(ICL.vec)
    vec.ICL <- c(vec.ICL,ICL.new)
    
    if (ICL.new > ICL.c) {
      w = which.max(ICL.vec)
      estim.new <- all_estim[[w]]
      estim.new$param_estim$Z <- lapply(1:data$Q,function(q){vq <- max.col(estim.new$param_estim$tau[[q]]); names(vq) <- data$names_ind[[q]]; return(vq)})
      
      RES[[niter.search + 1]] <- estim.new;  
      param.new <- estim.new$param_estim
      classif.new <- param.new$Z
      if (verbose) {
        mess <- paste(round(c(calc_vK(classif.new))),collapse = " " )
        print(paste("ICL :",round(ICL.new,2),". Nb of clusters: ", mess,sep = " "))
      }
    }
    
  }
  
  return(RES) 
}  

###########################"   




#####################  PLOT
plot_res = function(data,param_estim,mycol=NULL,thres=0.01){
  
  Q <-  data$Q
  vK_estim <- param_estim$vK
  type_node <- lapply(1:Q,function(q){rep(data$namesfg[q],vK_estim[q])})
  label_node <- lapply(1:Q,function(q){1:vK_estim[q]})
  if(is.null(mycol)){mycol <-  palette("default"); mycol <- mycol[-1]} 
  col_node <- lapply(1:Q,function(q){rep(mycol[q],vK_estim[q])})
  size_node <- lapply(1:Q,function(q){param_estim$lpi[[q]]})
  cum_vK =c(0,cumsum(vK_estim))
  code_node <- lapply(2:(Q+1),function(q){seq(cum_vK[q-1]+1,cum_vK[q],1)})
  
  
  N = sum(vK_estim)
  list_edges <- lapply(1:nrow(data$E),function(i){
    q.row = data$E[i,1]
    q.col = data$E[i,2]
    ltheta_i <- param_estim$ltheta[[i]]
    c1 <- rep(code_node[[q.row]],times=vK_estim[q.col])
    c2 <- rep(code_node[[q.col]],each = vK_estim[q.row])
    edges_i <- cbind(c1,c2,c(ltheta_i))
    edges_i <- as.data.frame(edges_i)
    edges_i$type <- rep(data$type_inter[i],length(c1))
    return(edges_i)})
  all_edges <- do.call("rbind", list_edges)
  
  # all_edges <- rbind(all_edges,c(2,2,10))
  # all_edges$type[nrow(all_edges)] = "diradj"
  w <- which(all_edges[,3]>thres)
  edges <- all_edges[w,c(1,2)]  
  curved <- 0.3*(all_edges[w,4]=="diradj") 
  curved <- runif(length(all_edges[w,4] == "diradj"),0,3)*(all_edges[w,4]=="diradj") 
  
  
  

  if(sum(data$type_inter=="diradj")==0){
    G<- make_empty_graph() + vertices(unlist(code_node))  
    G <- G %>%set_vertex_attr("label",value=unlist(label_node)) %>%set_vertex_attr("color",value=unlist(col_node))
    G <- G%>%set_vertex_attr("size",value=sqrt(unlist(size_node))*20+4)
    G <- G%>%set_edge_attr("width",value=all_edges[all_edges[,3]>=thres,3] )
    G <- G + edges(c(t(edges)))
    G <- as.undirected(G)
    
  }
  if(sum(data$type_inter=="diradj")>0){
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
  
   legend("left", c(data$namesfg), col=mycol[1:Q], border = "black", lty=1, lwd=4) 
}



plot_res_dattilo_LBM = function(data_dattilo, data_dattilo_all,param_estim_lbm_classique,mycol=NULL,thres=0.01){
  
  data <- data_dattilo_all
  Q <-  data$Q
  vK_estim <- param_estim_lbm_classique$vK
  type_node <- lapply(1:Q,function(q){rep(data$namesfg[q],vK_estim[q])})
  label_node <- lapply(1:Q,function(q){1:vK_estim[q]})
  if(is.null(mycol)){mycol <-  palette("default"); mycol <- mycol[-1]} 
  col_node <- lapply(1:Q,function(q){rep(mycol[q],vK_estim[q])})
  size_node <- lapply(1:Q,function(q){param_estim$lpi[[q]]})
  cum_vK =c(0,cumsum(vK_estim))
  code_node <- lapply(2:(Q+1),function(q){seq(cum_vK[q-1]+1,cum_vK[q],1)})
  
  
  
  
  N = sum(vK_estim)
  list_edges <- lapply(1:nrow(data$E),function(i){
    q.row = data$E[i,1]
    q.col = data$E[i,2]
    ltheta_i <- param_estim$ltheta[[i]]
    c1 <- rep(code_node[[q.row]],times=vK_estim[q.col])
    c2 <- rep(code_node[[q.col]],each = vK_estim[q.row])
    edges_i <- cbind(c1,c2,c(ltheta_i))
    edges_i <- as.data.frame(edges_i)
    edges_i$type <- rep(data$type_inter[i],length(c1))
    return(edges_i)})
  all_edges <- do.call("rbind", list_edges)
  w <- which(all_edges[,3]>thres)
  edges <- all_edges[w,c(1,2)]  
  curved <- 0.3*(all_edges[w,4]=="diradj") 
  
  
  
  if(sum(data$type_inter=="diradj")==0){
    G<- make_empty_graph() + vertices(unlist(code_node))  
    G <- G %>%set_vertex_attr("label",value=unlist(label_node)) %>%set_vertex_attr("color",value=unlist(col_node))
    G <- G%>%set_vertex_attr("size",value=sqrt(unlist(size_node))*20+4)
    G <- G%>%set_edge_attr("width",value=all_edges[all_edges[,3]>=thres,3] )
    G <- G + edges(c(t(edges)))
    G <- as.undirected(G)
    
  }
  if(sum(data$type_inter=="diradj")>0){
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
  
  legend("left", c(data$namesfg), col=mycol[1:Q], border = "black", lty=1, lwd=4) 
}



