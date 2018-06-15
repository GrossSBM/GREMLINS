estim_multi_init <- function(data,vK,nb_cores=NULL){
  #   
  
  os <- Sys.info()["sysname"]
  if ((os != 'Windows') & (is.null(nb_cores))) {nb_cores = detectCores(all.tests = FALSE, logical = TRUE) %/% 2}
  
  Q <- data$Q
  vdistrib <- data$vdistrib
  

  
  #------------------------  Initialisation of the number of clusters and the classification. 
  
  param.init <- genBMfit$new(vK = vK,vdistrib = vdistrib)

  classif.CAH <- initialize(data,param.init,method = "CAH")$groups
  
  #----------------------   Initialisation of the algorithm
  
   
  estim.0 <- data$estime(classif.CAH);
  param.0 <- estim.0$param_estim
  classif.0 <- lapply(1:Q,function(q){max.col(param.0$tau[[q]])})
  ICL.0 <- estim.0$ICL
  

  
  #----------------------------------------------  ALGORITHM
  
  
  
  ### step 1 -> M(+1)
  F_forward_q <-  function(q){
    classif_forward_q <- split_classif(classif.0,q,data,100)
    res <- do.call(c, list(classif_forward_q))
    return(res)}
  list_classif_init_forward = list()
  for (q in 1:Q) { list_classif_init_forward <- do.call(c,list(list_classif_init_forward,F_forward_q(q)))}
  
  
  if (os=="Windows"){
    all_estim_forward <- lapply(list_classif_init_forward,function(init){
      estim.c.l <- data$estime(init)})
  }else{
    all_estim_forward <- mclapply(list_classif_init_forward,function(init){estim.c.l <- data$estime(init)},mc.cores = nb_cores)
  }
  
  ### step 2 -> M(-1)
  F_backward_q <-  function(q){
    classif_backward_q <- merge_classif(classif.0,q,1)
    res <- do.call(c, list(classif_backward_q))
    return(res)}
  list_classif_init_backward = list()
  for (q in 1:Q) { list_classif_init_backward <- do.call(c,list(list_classif_init_backward,F_backward_q(q)))}
  
  
  if (os == "Windows") {
    all_estim_backward <- lapply(list_classif_init_backward,function(init){
    estim.c.l <- data$estime(init)})
  }else{
    all_estim_backward <- mclapply(list_classif_init_backward,function(init){estim.c.l <- data$estime(init)},mc.cores = nb_cores)
  }
  
  ICL.vec.forward <- vapply(all_estim_forward,function(u){u$ICL},1)
  ICL.vec.backward <- vapply(all_estim_backward,function(u){u$ICL},1)
  
  w.forward = which.max(ICL.vec.forward)
  estim.new.forward <- all_estim_forward[[w.forward]]
  param.new.forward <- estim.new.forward$param_estim
  classif.new.forward <- lapply(1:Q,function(q){max.col(param.new.forward$tau[[q]])})
  
  q_forward <- which(param.new.forward$vK != vK)
  init_forward <- merge_classif(classif.new.forward,q_forward,1)
  if (os == "Windows") {
    last_estim_forward <- lapply(init_forward,function(init){estim.c.l <- data$estime(init)})
  }else{
    last_estim_forward <- mclapply(init_forward,function(init){estim.c.l <- data$estime(init)},mc.cores = nb_cores)
  }
  
  
  ###########################"" 
  w.backward = which.max(ICL.vec.backward)
  estim.new.backward <- all_estim_backward[[w.backward]]
  param.new.backward <- estim.new.backward$param_estim
  classif.new.backward <- lapply(1:Q,function(q){max.col(param.new.backward$tau[[q]])})
  
  q_backward <- which(param.new.backward$vK != vK)
  init_backward <- split_classif(classif.new.backward,q_backward,data,100)
  if(os=="Windows"){
    last_estim_backward <- lapply(init_backward,function(init){estim.c.l <- data$estime(init)})
  }else{
    last_estim_backward <- mclapply(init_backward,function(init){estim.c.l <- data$estime(init)},mc.cores = nb_cores)
  }

  
  ICL.vec.forward <- vapply(last_estim_forward,function(u){u$ICL},1)
  ICL.vec.backward <- vapply(last_estim_backward,function(u){u$ICL},1)
  

  ### now we have find the best of the bests!!!!! 
  w.star.1 = which.max(ICL.vec.forward)
  w.star.2 = which.max(ICL.vec.backward)
  
  best = as.character(which.max(c(ICL.0,ICL.vec.forward[w.star.1],ICL.vec.backward[w.star.2])))
  
  res <- switch(best,
         "1" = {estim.0},
         "2" = {last_estim_forward[[w.star.1]]},
         "3" = {last_estim_backward[[w.star.2]]})
  res$param_estim$Z <- lapply(1:data$Q,function(q){vq <- max.col(res$param_estim$tau[[q]]); names(vq) <- data$names_ind[[q]]; return(vq)})
  res$classif <- lapply(res$param_estim$Z,function(z){sort(z)})
  return(res)
}  

###########################"   


