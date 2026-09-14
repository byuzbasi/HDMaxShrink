test_that("compiled solver converges without a Gram inverse", {
  set.seed(10)
  n <- 30L
  p <- 50L
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(1, -0.8, 0.6, rep(0, p - 3L))
  y <- drop(X %*% beta + rnorm(n, sd = 0.5))

  fit <- sqrt_lasso_hd(
    X,
    y,
    max_iter = 10000L,
    tol = 1e-4,
    check_every = 25L
  )
  expect_true(fit$converged)
  expect_length(fit$beta, p)
  expect_true(all(is.finite(fit$beta)))
  expect_lte(fit$kkt_residual, 1e-4)
})

test_that("zero-block projection is exact when q exceeds n", {
  set.seed(11)
  n <- 28L
  p <- 60L
  p1 <- 5L
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(0.8, p1), rep(0, p - p1))
  y <- drop(X %*% beta + rnorm(n))

  fit <- fit_hd_shrinkage(
    X,
    y,
    restricted_set = (p1 + 1L):p,
    test = "max",
    precision = "diagonal",
    assume_independent = TRUE,
    bootstrap_B = 99L,
    bootstrap_seed = 19L,
    max_iter = 10000L,
    tol = 1e-4
  )
  expect_gt(fit$dimensions[["q"]], fit$dimensions[["n"]])
  expect_equal(fit$beta$restricted[(p1 + 1L):p], rep(0, p - p1))
  expect_equal(
    fit$beta$restricted[seq_len(p1)],
    fit$beta$full[seq_len(p1)],
    tolerance = 1e-12
  )
  expect_equal(fit$restricted_solver$restriction_violation, 0)
  expect_identical(fit$restricted_solver$construction, "projection_of_full")
  expect_identical(fit$score_representation, "diagonal")
})

test_that("a general linear submodel is projected from the full fit", {
  set.seed(12)
  n <- 35L
  p <- 50L
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  C <- matrix(c(1, -1), nrow = 1L)
  target <- 0.25

  fit <- fit_hd_shrinkage(
    X,
    y,
    restricted_set = 1:2,
    C = C,
    t = target,
    endpoint = "regularized",
    test = "wald",
    precision = "diagonal",
    assume_independent = TRUE,
    max_iter = 10000L,
    tol = 1e-4
  )

  expect_equal(
    as.numeric(C %*% fit$beta$submodel[1:2]),
    target,
    tolerance = 1e-10
  )
  expect_equal(
    fit$beta$submodel[-(1:2)],
    fit$beta$full[-(1:2)],
    tolerance = 1e-12
  )
  expect_lte(fit$restricted_solver$restriction_violation, 1e-10)
})
