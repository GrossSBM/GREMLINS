
n_FG <- 3  #number of functional groups
v_K <- c(3,2,2) #number of clusters in each functional group
list_pi <- vector("list", 3); # parameters of clusterning
list_pi[[1]] <- c(0.4,0.3,0.3); list_pi[[2]] <- c(0.6,0.4); list_pi[[3]]  <- c(0.6,0.4)
E  = rbind(c(1,2),c(2,3),c(2,2),c(3,3))
v_distrib <- c('bernoulli','poisson','laplace','gaussian')
typeInter <- c( "inc", "inc"  ,  "adj" ,"diradj")
list_theta <- list()
list_theta[[1]] <- matrix(rbeta(v_K[E[1,1]] * v_K[E[1,2]],1.5,1.5 ),nrow = v_K[E[1,1]], ncol = v_K[E[1,2]] )
list_theta[[2]] <- matrix(rgamma(v_K[E[2,1]] * v_K[E[2,2]],7.5,1 ),nrow = v_K[E[2,1]], ncol = v_K[E[2,2]] )
list_theta[[3]] <- list()
list_theta[[3]] <- matrix(rnorm(v_K[E[3,1]] * v_K[E[3,2]],7.5,1 ),nrow = v_K[E[3,1]], ncol = v_K[E[3,2]] )
list_theta[[3]] <- 0.5*(list_theta[[3]] + t(list_theta[[3]])) # symetrisation for network 3
list_theta[[4]] <- list()
list_theta[[4]]$mean <- matrix(rnorm(v_K[E[4,1]] * v_K[E[4,2]],7.5,1 ),nrow = v_K[E[4,1]], ncol = v_K[E[4,2]] )
list_theta[[4]]$sd <- matrix(rgamma(v_K[E[4,1]] * v_K[E[4,2]],7.5,1 ),nrow = v_K[E[4,1]], ncol = v_K[E[4,2]] )
v_NQ = c(100,50,40)
dataSim <- rMBM(v_NQ ,E , typeInter, v_distrib, list_pi, list_theta, seed = NULL, namesFG= c('A','B','C'),keepClassif  = TRUE)
list_Net <- dataSim$list_Net

dataSim2 <- rMBM(v_NQ ,E , typeInter, v_distrib, list_pi, list_theta, seed = NULL, namesFG= c('A','B','C'),keepClassif  = FALSE)



test_that("The simulated object of correct type and  dimensions", {

  ########### dataSim must be a list
  expect_type(dataSim,"list")


  ########### list_Net must be a list
  expect_type(list_Net,"list")

  expect_equal(length(dataSim$list_Net),4)

  ######### Check the size of each matrix
  test_dim <- 1;
  for (net in 1:4) {
    dim_net <- dim(list_Net[[net]]$mat)
    dim_theo <- v_NQ[c(E[net,1],E[net,2])]
    test_dim <- test_dim * prod(dim_net == dim_theo)
  }
  expect_equal(test_dim,1)

  ######### Check the symmetry
  test_symmetry <- 1;
  w <- which(typeInter != 'inc')
  if (length(w) > 0) {
    for (net in w) {
      test_symmetry <- test_symmetry * (isSymmetric(list_Net[[net]]$mat) == (typeInter[net]=='adj'))
    }
  }
  expect_equal(test_symmetry,1)

  ######### Check the type of the content of each matrix
  test_type <- 1;
  for (net in 1:4) {

    if (v_distrib[net] == "bernoulli") {test_type <- test_type * all(list_Net[[net]]$mat %in% c(0,1)) }
    if (v_distrib[net] == "poisson") {test_type <- test_type * all(is.integer(list_Net[[net]]$mat)) }
    }
  expect_equal(test_type,1)
  }
)

test_that("We kept the classif or not",{
  expect_null(dataSim2$classif)
  expect_equal(length(dataSim$classif),n_FG)
}
)
