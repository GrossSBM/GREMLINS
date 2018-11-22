F_foreach_dopart <-
function(){
  cl <- makeCluster(3)
  registerDoParallel(cl)
  U <-  foreach(i = 1:100) %dopar% Func_test(i*10)
  stopCluster(cl)
  return(U)
}
