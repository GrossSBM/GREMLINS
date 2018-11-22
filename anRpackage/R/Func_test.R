Func_test <-
function(n){
  p = 6;
  X = matrix(rnorm(n*6)+3,n,p);
  res <- solve(t(X) %*% X)
  return(res)
}
