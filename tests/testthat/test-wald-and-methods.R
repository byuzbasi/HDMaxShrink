test_that("small-q Wald branch returns all classical estimators", {
  set.seed(30)
  n <- 35L
  p <- 55L
  X <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p)
  beta[6:9] <- c(1, -0.8, 0.7, -0.6)
  y <- drop(X %*% beta + rnorm(n))

  fit <- fit_hd_shrinkage(
    X,
    y,
    restricted_set = 1:3,
    test = "wald",
    precision = "diagonal",
    assume_independent = TRUE,
    max_iter = 10000L,
    tol = 1e-4
  )
  expect_identical(fit$test, "wald")
  expect_true(is.finite(fit$inference$W))
  expect_length(fit$beta$stein, p)
  expect_length(fit$beta$positive_part, p)
  expect_length(coef(fit, method = "PS"), p + 1L)
  expect_length(predict(fit, X[1:4, , drop = FALSE]), 4L)
})

test_that("diagonal precision requires an explicit independence assumption", {
  X <- matrix(rnorm(20 * 30), 20, 30)
  expect_error(
    estimate_precision_rows(X, 1:3, method = "diagonal"),
    "assume_independent"
  )
})
