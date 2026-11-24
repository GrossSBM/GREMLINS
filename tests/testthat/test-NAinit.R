context("Initialisation with missing values")

simuNA <- function(propNA = 0.2) {
  E <- rbind(c(1, 2), c(2, 2))
  sim <- rMBM(v_NQ = c(40, 40), E, c("inc", "diradj"), c("poisson", "bernoulli"),
              list(c(0.5, 0.5), c(0.3, 0.7)),
              list(matrix(c(6.1, 8.9, 6.6, 3), 2, 2), matrix(c(0.7, 1.0, 0.4, 0.6), 2, 2)),
              namesFG = c("A", "B"), seed = 2, keepClassif = TRUE)
  set.seed(1)
  for (e in 1:2) {
    m <- sim$list_Net[[e]]$mat
    m[sample(length(m), propNA * length(m))] <- NA
    # two individuals without any common observed entry
    m[1, 1:20] <- NA; m[2, 21:40] <- NA
    sim$list_Net[[e]]$mat <- m
  }
  sim
}

test_that("NA-safe distance never returns NA", {
  x <- matrix(c(NA, NA, 1, 0, 1, 1, NA, NA), 2)
  d <- GREMLINS:::distNA(x)
  expect_false(anyNA(d))
})

test_that("adjustClassif reaches the required number of blocks", {
  sim <- simuNA()
  dataR6 <- GREMLINS:::formattingData(sim$list_Net, c("poisson", "bernoulli"))
  z <- rep(1:4, length.out = 40)
  expect_equal(max(GREMLINS:::adjustClassif(z, 2, dataR6, 1)), 2)
  expect_equal(max(GREMLINS:::adjustClassif(z, 6, dataR6, 1)), 6)
})

test_that("multipartiteBM works with NA, with blockmodels and CAH initialisations", {
  sim <- simuNA()
  for (bm in c(TRUE, FALSE)) {
    res <- multipartiteBM(sim$list_Net, c("poisson", "bernoulli"), namesFG = c("A", "B"),
                          v_Kmin = 1, v_Kmax = 4, initBM = bm, verbose = FALSE)
    expect_equal(unname(res$fittedModel[[1]]$paramEstim$v_K), c(2, 2))
  }
})

test_that("multipartiteBMFixedModel works with NA and K_q = 1", {
  sim <- simuNA()
  for (vK in list(c(2, 2), c(1, 1), c(2, 1))) {
    res <- multipartiteBMFixedModel(sim$list_Net, c("poisson", "bernoulli"), namesFG = c("A", "B"),
                                    v_K = vK, verbose = FALSE)
    expect_true(all(res$fittedModel[[1]]$paramEstim$v_K <= vK))
  }
})

test_that("v_Kinit is reordered like v_Kmin / v_Kmax (3 FGs, circular permutation)", {
  set.seed(1)
  # internal order of the FGs is B, C, A
  net1 <- defineNetwork(matrix(rbinom(900, 1, .3), 30, 30), "inc", "B", "C")
  net2 <- defineNetwork(matrix(rbinom(900, 1, .3), 30, 30), "inc", "C", "A")
  # only the starting point is checked: few iterations, the final non-convergence error is ignored
  out <- capture.output(try(
    multipartiteBM(list(net1, net2), c("bernoulli", "bernoulli"), namesFG = c("A", "B", "C"),
                   v_Kinit = c(1, 2, 3), v_Kmax = 4, initBM = FALSE, verbose = TRUE, maxiterVEM = 2),
    silent = TRUE
  ))
  expect_true(any(grepl("starting from \\[ 2 3 1 \\]", out)))
})
