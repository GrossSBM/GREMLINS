context("testing clusters recovery")


set.seed(123)
v_K <- c(2,2)
n_FG <- 2
list_pi <- vector("list", 2);
list_pi[[1]] <- c(0.5,.5); list_pi[[2]] <- c(0.6,0.4)
E  = rbind(c(1,2),c(2,2))
v_distrib <- c('bernoulli','poisson')
typeInter <- c( "inc", "adj" )
list_theta <- list()
list_theta[[1]] <- matrix(rbeta(v_K[E[1,1]] * v_K[E[1,2]],1.5,1.5 ),nrow = v_K[E[1,1]], ncol = v_K[E[1,2]])
list_theta[[2]] <- matrix(rgamma(v_K[E[2,1]] * v_K[E[2,2]],7.5,1 ),nrow = v_K[E[2,1]], ncol = v_K[E[2,2]])
v_NQ = c(50,40)
dataSim <-  rMBM(v_NQ ,E , typeInter, v_distrib, list_pi, list_theta, seed=NULL, namesFG= c('A','B'),keepClassif = TRUE)
list_Net <- dataSim$list_Net
res <- multipartiteBM(list_Net, v_distrib = c("bernoulli","poisson"), namesFG = NULL, v_Kmin = 1,v_Kmax = 5,
                      v_Kinit = NULL,nbCores  = 2, verbose = TRUE, save=FALSE, maxiterVE = NULL)



library(aricode)
test_that("correct clustering 1st FG", {
  expect_gt(ARI(res$fittedModel[[1]]$paramEstim$Z[[1]],dataSim$classif[[1]]),.4)
})


test_that("correct clustering 2nd FG", {
  expect_gt(ARI(res$fittedModel[[1]]$paramEstim$Z[[2]],dataSim$classif[[2]]),.4)
})


clust <- extractClustersMBM (res,whichModel = 1)
test_that('extractClustersMBM',{
  expect_equal(length(clust),n_FG)
})


