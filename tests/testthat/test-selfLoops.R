context("Self-loops")

test_that("the diagonal is ignored by default and kept with selfLoops = TRUE", {
  set.seed(1)
  A <- matrix(rbinom(400, 1, .1), 20, 20)
  estim <- function(A, selfLoops = FALSE) {
    res <- multipartiteBMFixedModel(list(defineNetwork(A, "diradj", "FG", "FG", selfLoops = selfLoops)),
                                    v_distrib = "bernoulli", v_K = 1, verbose = FALSE)
    c(res$fittedModel[[1]]$paramEstim$list_theta[[1]])
  }
  A0 <- A; diag(A0) <- 0
  A1 <- A; diag(A1) <- 1
  expect_true(all(is.na(diag(defineNetwork(A, "diradj", "FG", "FG")$mat))))
  expect_equal(estim(A0), estim(A1))
  expect_equal(estim(A0), sum(A0) / 380, tolerance = 1e-6)
  expect_equal(estim(A1, selfLoops = TRUE), sum(A1) / 400, tolerance = 1e-6)
})

test_that("networks given without defineNetwork also ignore the diagonal", {
  set.seed(2)
  A <- matrix(rbinom(400, 1, .3), 20, 20); A <- 1 * ((A + t(A)) > 0)
  dataR6 <- GREMLINS:::formattingData(list(list(mat = A, typeInter = "adj", rowFG = "FG", colFG = "FG")),
                                      "bernoulli")
  expect_true(all(is.na(diag(dataR6$mats[[1]]))))
})

test_that("rMBM: no self-loops, incidence matrices untouched", {
  E <- rbind(c(1, 2), c(2, 2))
  sim <- rMBM(v_NQ = c(30, 30), E, c("inc", "diradj"), c("bernoulli", "poisson"),
              list(1, 1), list(matrix(1 - 1e-12), matrix(3)), seed = 1)
  expect_false(anyNA(sim$list_Net[[1]]$mat))
  expect_true(all(diag(sim$list_Net[[1]]$mat) == 1))
  expect_true(all(is.na(diag(sim$list_Net[[2]]$mat))))
})
