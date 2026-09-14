test_that("dual Ridge agrees with the direct n-dimensional formula", {
  set.seed(50)
  n <- 14L
  p <- 25L
  X <- matrix(rnorm(n * p), n, p)
  X <- sweep(X, 2L, colMeans(X), "-")
  X <- sweep(X, 2L, sqrt(colMeans(X^2)), "/")
  y <- rnorm(n)
  lambda <- c(0.15, 0.8)

  observed <- HDMaxShrink:::cpp_ridge_dual_path(X, y, lambda)
  for (index in seq_along(lambda)) {
    expected <- drop(t(X) %*% solve(
      X %*% t(X) + n * lambda[index] * diag(n), y
    ))
    expect_equal(
      observed$beta[, index], expected,
      tolerance = 1e-11
    )
  }
  expect_true(all(is.finite(observed$gcv)))
  expect_equal(
    observed$gcv,
    (observed$rss / n) /
      (1 - (observed$degrees_freedom + 1) / n)^2,
    tolerance = 1e-12
  )
})

test_that("the MCP kernel matches its standardized univariate update", {
  x <- seq_len(12L) - mean(seq_len(12L))
  x <- x / sqrt(mean(x^2))
  z <- 1.5
  lambda <- 0.5
  gamma <- 3
  expected <- (z - lambda) / (1 - 1 / gamma)

  fit <- HDMaxShrink:::cpp_mcp_path(
    X = matrix(x, ncol = 1L),
    y = x * z,
    lambda = lambda,
    beta_init = 0,
    gamma = gamma,
    max_iter = 100L,
    tol = 1e-12,
    zero_tol = 0
  )
  expect_true(fit$converged[1L])
  expect_equal(fit$beta[1L, 1L], expected, tolerance = 1e-11)
})

test_that("MCP uses all X and is not forced into the tested null", {
  set.seed(51)
  n <- 50L
  p <- 80L
  tested <- 4:p
  X <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p)
  beta[c(1L, 2L, 5L)] <- c(1.2, -1, 2)
  y <- drop(X %*% beta + rnorm(n, sd = 0.15))

  fit <- fit_ridge_mcp_shrinkage(
    X,
    y,
    tested_set = tested,
    ridge_lambda = 0.5,
    mcp_lambda = 0.08,
    precision = "diagonal",
    assume_independent = TRUE,
    bootstrap_B = 99L,
    bootstrap_seed = 52L,
    pilot_max_iter = 10000L,
    pilot_tol = 1e-4
  )

  expect_identical(fit$endpoint, "ridge_mcp")
  expect_false(fit$restriction$exact_endpoint_restriction)
  expect_identical(
    fit$restriction$construction,
    "data_adaptive_MCP_all_X"
  )
  expect_gt(abs(fit$beta$submodel[5L]), 0.5)
  expect_true(5L %in% fit$restricted_solver$support)
  expect_gt(fit$restricted_solver$restriction_violation, 0.5)
  expect_equal(fit$beta$full, fit$beta$ridge)
  expect_equal(fit$beta$submodel, fit$beta$mcp)
  expect_equal(unname(coef(fit, method = "RIDGE")[-1L]), fit$beta$full)
  expect_equal(unname(coef(fit, method = "MCP")[-1L]), fit$beta$submodel)
  expect_equal(
    fit$beta$positive_part,
    fit$beta$submodel + fit$shrinkage$positive_weight *
      (fit$beta$full - fit$beta$submodel),
    tolerance = 1e-12
  )
  expect_equal(
    fit$beta$protected_positive_part,
    if (fit$inference$reject) fit$beta$full else fit$beta$positive_part,
    tolerance = 1e-12
  )
})

test_that("automatic Ridge and MCP tuning return finite endpoints", {
  set.seed(53)
  n <- 35L
  p <- 60L
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(1.1, -0.9, 0.7, rep(0, p - 3L))
  y <- drop(X %*% beta + rnorm(n))

  fit <- fit_ridge_mcp_shrinkage(
    X,
    y,
    tested_set = 4:p,
    ridge_lambda_grid = 10^seq(-2, 1, length.out = 8L),
    mcp_nlambda = 15L,
    precision = "diagonal",
    assume_independent = TRUE,
    bootstrap_B = 99L,
    bootstrap_seed = 54L,
    pilot_max_iter = 10000L,
    pilot_tol = 1e-4
  )

  expect_identical(fit$full_solver$selection, "GCV")
  expect_identical(fit$restricted_solver$selection, "EBIC")
  expect_true(fit$restricted_solver$converged)
  expect_true(all(is.finite(fit$beta$full)))
  expect_true(all(is.finite(fit$beta$submodel)))
  expect_lte(length(fit$restricted_solver$support), p)
})
