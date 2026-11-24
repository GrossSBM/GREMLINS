context("Order of the functional groups")

# A has 1 block, B 2 blocks, C 3 blocks. Internal order of the FGs: B, C, A (circular permutation of namesFG)
simuOrder <- function() {
  set.seed(2)
  zA <- rep(1, 30); zB <- rep(1:2, each = 20); zC <- rep(1:3, each = 20)
  mk <- function(z1, z2, p) {
    m <- matrix(rbinom(length(z1) * length(z2), 1, p[cbind(rep(z1, length(z2)), rep(z2, each = length(z1)))]), length(z1))
    rownames(m) <- paste0("r", seq_along(z1)); colnames(m) <- paste0("c", seq_along(z2)); m
  }
  net1 <- defineNetwork(mk(zB, zC, matrix(c(.9, .1, .5, .1, .9, .1), 2, 3, byrow = TRUE)), "inc", "B", "C")
  net2 <- defineNetwork(mk(zC, zA, matrix(c(.8, .2, .5), 3, 1)), "inc", "C", "A")
  list(list_Net = list(net1, net2), z = list(A = zA, B = zB, C = zC))
}

expectGoodLabels <- function(fit) {
  p <- fit$fittedModel[[1]]$paramEstim
  expect_equal(as.numeric(p$v_K[c("A", "B", "C")]), c(1, 2, 3))
  expect_equal(names(p$Z), names(p$v_K))
  expect_equal(as.numeric(sapply(p$Z, max)), as.numeric(p$v_K))
}

test_that("outputs are correctly labelled", {
  s <- simuOrder()
  expectGoodLabels(multipartiteBM(s$list_Net, c("bernoulli", "bernoulli"), namesFG = c("A", "B", "C"), v_Kmax = 5, verbose = FALSE))
  expectGoodLabels(multipartiteBMFixedModel(s$list_Net, c("bernoulli", "bernoulli"), namesFG = c("A", "B", "C"), v_K = c(1, 2, 3), verbose = FALSE))
  expectGoodLabels(multipartiteBMFixedModel(s$list_Net, c("bernoulli", "bernoulli"), namesFG = c("A", "B", "C"), v_K = c(1, 2, 3), verbose = FALSE, initBM = FALSE))
})

test_that("givenclassif / classifInit follow namesFG or their names", {
  s <- simuOrder()
  gc <- unname(s$z)                       # order of namesFG
  gcNamed <- s$z[c("C", "A", "B")]        # named, any order
  for (g in list(gc, gcNamed)) {
    expectGoodLabels(multipartiteBM(s$list_Net, c("bernoulli", "bernoulli"), namesFG = c("A", "B", "C"), v_Kmax = 5, givenclassif = g, verbose = FALSE))
    expectGoodLabels(multipartiteBMFixedModel(s$list_Net, c("bernoulli", "bernoulli"), namesFG = c("A", "B", "C"), v_K = NULL, classifInit = g, verbose = FALSE))
  }
  expect_error(multipartiteBM(s$list_Net, c("bernoulli", "bernoulli"), namesFG = c("A", "B", "C"), givenclassif = list(s$z$B, s$z$A, s$z$C), verbose = FALSE), "must be of length")
  expect_error(multipartiteBM(s$list_Net, c("bernoulli", "bernoulli"), namesFG = c("A", "B", "C"), givenclassif = list(X = 1, Y = 1, Z = 1), verbose = FALSE), "do not match")
})
