test_that("all-X FM is fitted once and is independent of the A/B partition", {
  set.seed(20260824)
  n <- 48L
  p <- 72L
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(1.1, -0.9, 0.7, 0.6, rep(0, p - 4L))
  y <- drop(X %*% beta + rnorm(n, sd = 0.8))
  gaussian_draws <- matrix(rnorm(n * 99L), n, 99L)

  first <- fit_partial_sqrt_lasso_shrinkage(
    X = X,
    y = y,
    core_set = 1:3,
    tested_set = 4:p,
    lambda = 0.16,
    bootstrap_B = 99L,
    gaussian_draws = gaussian_draws,
    max_iter = 10000L,
    tol = 1e-6,
    full_endpoint = "all_x_square_root_lasso"
  )
  second <- fit_partial_sqrt_lasso_shrinkage(
    X = X,
    y = y,
    core_set = 1:5,
    tested_set = 6:p,
    lambda = 0.16,
    bootstrap_B = 99L,
    gaussian_draws = gaussian_draws,
    max_iter = 10000L,
    tol = 1e-6,
    full_endpoint = "all_x_square_root_lasso"
  )

  expect_identical(first$full_endpoint, "all_x_square_root_lasso")
  expect_identical(first$full_solver$construction, "all_x_square_root_lasso")
  expect_identical(first$full_solver$penalty_scope, "all_X")
  expect_true(first$full_solver$standardized_response)
  expect_true(all(unlist(first$standardized_intercept) == 0))
  expect_equal(first$beta$full, second$beta$full, tolerance = 1e-10)
  expect_equal(first$intercept$full, second$intercept$full, tolerance = 1e-10)
  expect_equal(first$beta$submodel[4:p], numeric(p - 3L), tolerance = 0)
  expect_equal(second$beta$submodel[6:p], numeric(p - 5L), tolerance = 0)
  expect_equal(length(first$full_solver$penalty_loadings), p)
  expect_true(all(first$full_solver$penalty_loadings == 1))

  prep <- HDMaxShrink:::.prepare_design(X, y, scale_y = TRUE)
  direct <- HDMaxShrink:::.sqrt_lasso_standardized(
    X = prep$X,
    y = prep$y,
    lambda = 0.16,
    penalty_factor = rep(1, p),
    max_iter = 10000L,
    tol = 1e-6
  )
  direct_beta <- as.numeric(direct$beta) * prep$y_scale / prep$x_scale
  direct_intercept <- prep$y_center - sum(prep$x_center * direct_beta)
  expect_equal(first$beta$full, direct_beta, tolerance = 1e-10)
  expect_equal(first$intercept$full, direct_intercept, tolerance = 1e-10)
})

test_that("all-X endpoint has an unambiguous loading and scale contract", {
  set.seed(20260825)
  X <- matrix(rnorm(30 * 45), 30, 45)
  y <- rnorm(30)

  expect_error(
    fit_partial_sqrt_lasso_shrinkage(
      X, y, core_set = 1:3,
      penalty_loadings = rep(1, 42),
      bootstrap_B = 99,
      full_endpoint = "all_x_square_root_lasso"
    ),
    "all_x_penalty_loadings"
  )
  expect_error(
    fit_partial_sqrt_lasso_shrinkage(
      X, y, core_set = 1:3,
      all_x_penalty_loadings = rep(1, 44),
      bootstrap_B = 99,
      full_endpoint = "all_x_square_root_lasso"
    ),
    "p finite positive"
  )
  expect_error(
    fit_partial_sqrt_lasso_shrinkage(
      X, y, core_set = 1:3,
      standardize_y = FALSE,
      bootstrap_B = 99,
      full_endpoint = "all_x_square_root_lasso"
    ),
    "standardize_y must be TRUE"
  )
})
