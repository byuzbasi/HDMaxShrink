test_that("profiled square-root LASSO is a correlated-design FM endpoint", {
  set.seed(81)
  n <- 55L
  p <- 90L
  core <- 1:5
  tested <- 6:p
  common <- matrix(rnorm(n * 3L), n, 3L)
  loadings <- matrix(rnorm(3L * p), 3L, p)
  X <- sqrt(0.65) * common %*% loadings / sqrt(3) +
    sqrt(0.35) * matrix(rnorm(n * p), n, p)
  beta <- numeric(p)
  beta[c(1:5, 8, 13)] <- c(1.2, -1, 0.8, 0.6, -0.5, 1.4, -1.1)
  y <- drop(X %*% beta + rnorm(n, sd = 0.7))
  gaussian_draws <- matrix(rnorm(n * 99L), n, 99L)

  fit <- fit_partial_sqrt_lasso_shrinkage(
    X,
    y,
    core_set = core,
    tested_set = tested,
    full_endpoint = "profiled_square_root_lasso",
    gaussian_draws = gaussian_draws,
    max_iter = 20000L,
    tol = 1e-5
  )

  expect_identical(fit$endpoint, "partial_sqrt_lasso_null")
  expect_identical(fit$full_endpoint, "profiled_square_root_lasso")
  expect_identical(
    fit$full_solver$construction,
    "profiled_partial_square_root_lasso"
  )
  expect_identical(fit$test, "max_partial_t")
  expect_identical(fit$test_calibration, "conditional_gaussian")
  expect_identical(fit$inference$type, "gaussian_max_partial_t")
  expect_identical(fit$threshold$type, "not_applicable")
  expect_identical(fit$threshold$studentization, "not_applied")
  expect_false(fit$precision$required)
  expect_false(fit$precision$assume_independent)
  expect_equal(fit$beta$full, fit$beta$full_regularized, tolerance = 0)
  expect_equal(
    fit$beta$full[tested],
    as.numeric(fit$full_solver$beta_profiled) /
      fit$preprocessing$x_scale[tested],
    tolerance = 1e-12
  )
  expect_equal(fit$beta$submodel[tested], numeric(length(tested)))
  expect_identical(fit$restriction$construction, "exact_null_refit_svd")
  expect_true(fit$restriction$exact_endpoint_restriction)
  expect_null(fit$beta$full_debiased)
  expect_null(fit$inference$theta_tilde)
  expect_error(coef(fit, method = "FDSL"), "undefined")
  expect_equal(
    unname(coef(fit, method = "FM")[-1L]),
    unname(coef(fit, method = "FSL")[-1L]),
    tolerance = 0
  )
  expect_equal(
    fit$beta$positive_part,
    fit$beta$submodel + fit$shrinkage$positive_weight *
      (fit$beta$full - fit$beta$submodel),
    tolerance = 1e-12
  )
  expect_equal(
    as.numeric(predict(fit, X[1:4, , drop = FALSE], method = "PS")),
    as.numeric(
      fit$intercept$positive_part +
        X[1:4, , drop = FALSE] %*% fit$beta$positive_part
    ),
    tolerance = 1e-12
  )

  printed <- capture.output(print(fit))
  expect_true(any(grepl("profiled partial square-root-LASSO", printed)))
  expect_true(any(grepl("active tested=", printed)))
})

test_that("profiled endpoint rejects incompatible calibration arguments", {
  set.seed(82)
  X <- matrix(rnorm(35 * 60), 35, 60)
  y <- rnorm(35)

  expect_error(
    fit_partial_sqrt_lasso_shrinkage(
      X,
      y,
      core_set = 1:4,
      full_endpoint = "profiled_square_root_lasso",
      test_calibration = "residual_multiplier"
    ),
    "requires.*conditional_gaussian"
  )
  expect_error(
    fit_partial_sqrt_lasso_shrinkage(
      X,
      y,
      core_set = 1:4,
      threshold_eta = 1,
      full_endpoint = "profiled_square_root_lasso"
    ),
    "threshold_eta is not used"
  )
})

test_that("max-score wrapper exposes the profiled FM without changing defaults", {
  set.seed(83)
  n <- 40L
  p <- 65L
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  gaussian_draws <- matrix(rnorm(n * 99L), n, 99L)

  profiled <- max_score_test_hd(
    X,
    y,
    core_set = 1:4,
    gaussian_draws = gaussian_draws,
    full_endpoint = "profiled_square_root_lasso",
    max_iter = 10000L,
    tol = 1e-4
  )
  default <- fit_partial_sqrt_lasso_shrinkage(
    X,
    y,
    core_set = 1:4,
    assume_independent = TRUE,
    gaussian_draws = gaussian_draws,
    max_iter = 10000L,
    tol = 1e-4
  )

  expect_identical(profiled$full_endpoint, "profiled_square_root_lasso")
  expect_identical(
    profiled$full_solver$construction,
    "profiled_partial_square_root_lasso"
  )
  expect_equal(
    profiled$full_model_beta,
    profiled$full_solver$beta / default$preprocessing$x_scale,
    tolerance = 1e-12
  )
  expect_null(default$full_endpoint)
  expect_identical(
    default$full_solver$construction,
    "thresholded_partial_debiased_square_root_lasso"
  )
  expect_identical(default$threshold$type, "marginal_gaussian_pfer")
})
