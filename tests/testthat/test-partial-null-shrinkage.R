test_that("SVD null geometry reproduces a QR core refit", {
  set.seed(70)
  n <- 45L
  p1 <- 5L
  q <- 23L
  X_core <- matrix(rnorm(n * p1), n, p1)
  X_tested <- matrix(rnorm(n * q), n, q)
  y <- rnorm(n)

  geometry <- HDMaxShrink:::cpp_partial_null_geometry(X_core, X_tested)
  state <- HDMaxShrink:::cpp_partial_null_apply(
    y,
    geometry$coefficient_operator,
    geometry$orthonormal_basis,
    geometry$residualized_tested
  )

  expected <- qr.coef(qr(X_core), y)
  expect_equal(as.numeric(state$beta_core), as.numeric(expected), tolerance = 1e-11)
  expect_equal(
    crossprod(X_core, geometry$residualized_tested),
    matrix(0, p1, q),
    tolerance = 1e-10
  )
  expect_equal(geometry$rank, p1)
  expect_lt(geometry$orthogonality_error, 1e-10)
})

test_that("partial one-step scores equal the direct formulas", {
  set.seed(71)
  n <- 35L
  q <- 40L
  V <- matrix(rnorm(n * q), n, q)
  y <- rnorm(n)
  beta <- rnorm(q, sd = 0.1)

  observed <- HDMaxShrink:::cpp_partial_debiased_scores(y, V, beta)
  residual <- drop(y - V %*% beta)
  second_moment <- colMeans(V^2)
  correction <- drop(crossprod(V, residual)) / n / second_moment
  raw_psi <- sweep(V * residual, 2L, second_moment, "/")
  centered_psi <- sweep(raw_psi, 2L, colMeans(raw_psi), "-")

  expect_equal(
    as.numeric(observed$theta_tilde), beta + correction, tolerance = 1e-12
  )
  expect_equal(observed$psi_centered, centered_psi, tolerance = 1e-12)
  expect_equal(
    as.numeric(observed$variance), colMeans(centered_psi^2), tolerance = 1e-12
  )
})

test_that("FM is thresholded partial-debiased and SM is an exact-null refit", {
  set.seed(72)
  n <- 60L
  p <- 120L
  core <- 1:4
  tested <- 5:p
  X <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p)
  beta[1:5] <- c(1.2, -1, 0.8, 0.6, 4)
  y <- drop(X %*% beta + rnorm(n, sd = 0.5))
  multipliers <- matrix(rnorm(n * 99L), n, 99L)

  fit <- fit_partial_sqrt_lasso_shrinkage(
    X,
    y,
    core_set = core,
    tested_set = tested,
    lambda = 0.15,
    assume_independent = TRUE,
    bootstrap_multipliers = multipliers,
    max_iter = 20000L,
    tol = 1e-5
  )

  expect_identical(fit$endpoint, "partial_sqrt_lasso_null")
  expect_identical(fit$test, "max_partial_t")
  expect_identical(fit$test_calibration, "conditional_gaussian")
  expect_identical(fit$inference$type, "gaussian_max_partial_t")
  expect_equal(
    fit$inference$statistic,
    fit$inference$max_partial_t,
    tolerance = 0
  )
  expect_true(is.finite(fit$inference$score_statistic))
  expect_identical(
    fit$inference$shrinkage_calibration_type,
    "inverse_moment"
  )
  expect_equal(
    fit$inference$null_energy,
    fit$inference$null_second_moment,
    tolerance = 0
  )
  expect_equal(
    fit$shrinkage$shrinkage_calibration,
    fit$inference$null_inverse_moment_calibration,
    tolerance = 0
  )
  expect_identical(fit$restriction$construction, "exact_null_refit_svd")
  expect_true(fit$restriction$exact_endpoint_restriction)
  expect_equal(fit$beta$submodel[tested], numeric(length(tested)))
  residual_df <- n - fit$restricted_solver$rank - 1L
  z_score <- sqrt(residual_df) * fit$inference$theta_tilde /
    fit$inference$standard_error
  retained <- abs(z_score) > fit$threshold$value
  expect_equal(fit$threshold$eta, 5)
  expect_equal(fit$threshold$residual_df, residual_df)
  expect_equal(
    fit$threshold$studentization_correction,
    sqrt(residual_df / n),
    tolerance = 0
  )
  expect_identical(
    fit$threshold$studentization,
    "residual_df_adjusted_partial_debiased"
  )
  expect_equal(
    fit$threshold$value,
    qnorm(1 - fit$threshold$eta / (2 * length(tested))),
    tolerance = 1e-14
  )
  expect_identical(fit$threshold$retained, retained)
  expect_equal(
    fit$beta$full[tested],
    fit$inference$theta_tilde * retained,
    tolerance = 1e-12
  )
  expect_equal(
    fit$beta$full_regularized[tested],
    as.numeric(fit$full_solver$beta_profiled) /
      fit$preprocessing$x_scale[tested],
    tolerance = 1e-12
  )
  expect_gt(abs(fit$beta$full[5L]), 3)
  expect_gt(max(abs(fit$beta$full[core] - fit$beta$submodel[core])), 0.1)
  expect_gt(fit$inference$statistic, fit$inference$critical_value)
  expect_equal(
    fit$beta$positive_part,
    fit$beta$submodel + fit$shrinkage$positive_weight *
      (fit$beta$full - fit$beta$submodel),
    tolerance = 1e-12
  )
  expect_equal(
    fit$beta$protected_positive_part,
    if (fit$inference$reject) fit$beta$full else fit$beta$positive_part,
    tolerance = 0
  )
  expect_equal(
    fit$shrinkage$protected_weight,
    if (fit$inference$reject) 1 else fit$shrinkage$positive_weight,
    tolerance = 0
  )
  expect_equal(
    fit$beta$preliminary_test,
    if (fit$inference$reject) fit$beta$full else fit$beta$submodel,
    tolerance = 0
  )
  expect_equal(
    unname(coef(fit, method = "SM")[-1L]),
    fit$beta$submodel
  )
  expect_equal(
    unname(coef(fit, method = "PPS")[-1L]),
    fit$beta$protected_positive_part
  )
  expect_equal(
    as.numeric(predict(fit, X[1:3, , drop = FALSE], method = "PPS")),
    as.numeric(
      fit$intercept$protected_positive_part +
        X[1:3, , drop = FALSE] %*% fit$beta$protected_positive_part
    ),
    tolerance = 1e-12
  )
})

test_that("threshold studentization uses residual degrees of freedom", {
  set.seed(77)
  n <- 38L
  p <- 65L
  core <- 1:6
  tested <- 7:p
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  fit <- fit_partial_sqrt_lasso_shrinkage(
    X,
    y,
    core_set = core,
    tested_set = tested,
    threshold_eta = 0.5,
    assume_independent = TRUE,
    bootstrap_B = 99L,
    bootstrap_seed = 77L,
    max_iter = 10000L,
    tol = 1e-4
  )

  residual_df <- n - fit$null_geometry$rank - 1L
  expected <- sqrt(residual_df) * fit$inference$theta_tilde /
    fit$inference$standard_error
  old_scale <- sqrt(n) * fit$inference$theta_tilde /
    fit$inference$standard_error

  expect_equal(fit$threshold$studentized, expected, tolerance = 1e-13)
  expect_equal(
    fit$threshold$studentized,
    sqrt(residual_df / n) * old_scale,
    tolerance = 1e-13
  )
  expect_equal(fit$full_solver$threshold_residual_df, residual_df)
  expect_equal(
    fit$full_solver$threshold_studentization_correction,
    sqrt(residual_df / n),
    tolerance = 0
  )
})

test_that("fixed primary and vanishing sensitivity eta rules are distinct", {
  expect_equal(HDMaxShrink:::.default_threshold_eta(980L), 5)
  expect_equal(
    HDMaxShrink:::.vanishing_threshold_eta(980L),
    1 / log(980),
    tolerance = 0
  )
  expect_equal(HDMaxShrink:::.default_threshold_eta(1L), 0.5)
  expect_equal(HDMaxShrink:::.vanishing_threshold_eta(1L), 0.5)
})

test_that("legacy residual-multiplier and second-moment mode is explicit", {
  set.seed(76)
  n <- 42L
  p <- 75L
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  multipliers <- matrix(rnorm(n * 99L), n, 99L)

  fit <- fit_partial_sqrt_lasso_shrinkage(
    X,
    y,
    core_set = 1:5,
    assume_independent = TRUE,
    bootstrap_multipliers = multipliers,
    test_calibration = "residual_multiplier",
    shrinkage_calibration = "second_moment",
    max_iter = 10000L,
    tol = 1e-4
  )

  expect_identical(fit$test, "max")
  expect_identical(fit$test_calibration, "residual_multiplier")
  expect_identical(fit$inference$type, "max")
  expect_true(is.na(fit$inference$max_partial_t))
  expect_equal(
    fit$inference$shrinkage_calibration,
    fit$inference$null_second_moment,
    tolerance = 0
  )
  expect_identical(
    fit$shrinkage$shrinkage_calibration_type,
    "second_moment"
  )
  expect_equal(
    fit$beta$positive_part,
    fit$beta$submodel + fit$shrinkage$positive_weight *
      (fit$beta$full - fit$beta$submodel),
    tolerance = 1e-12
  )
})

test_that("cached geometry is accepted only for the same design", {
  set.seed(73)
  n <- 40L
  p <- 70L
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  first <- fit_partial_sqrt_lasso_shrinkage(
    X,
    y,
    core_set = 1:5,
    assume_independent = TRUE,
    bootstrap_B = 99L,
    bootstrap_seed = 1L,
    max_iter = 10000L,
    tol = 1e-4
  )
  second <- fit_partial_sqrt_lasso_shrinkage(
    X,
    y,
    core_set = 1:5,
    assume_independent = TRUE,
    bootstrap_B = 99L,
    bootstrap_seed = 1L,
    max_iter = 10000L,
    tol = 1e-4,
    null_geometry = first$null_geometry,
    full_operator_norm = first$full_solver$operator_norm
  )
  expect_equal(first$beta$full, second$beta$full, tolerance = 1e-12)

  changed <- X
  changed[1L, 1L] <- changed[1L, 1L] + 0.25
  expect_error(
    fit_partial_sqrt_lasso_shrinkage(
      changed,
      y,
      core_set = 1:5,
      assume_independent = TRUE,
      bootstrap_B = 99L,
      null_geometry = first$null_geometry
    ),
    "different design"
  )
})

test_that("partial framework enforces its declared design assumptions", {
  set.seed(74)
  X <- matrix(rnorm(30 * 50), 30, 50)
  y <- rnorm(30)
  expect_error(
    fit_partial_sqrt_lasso_shrinkage(X, y, core_set = 1:3),
    "assume_independent"
  )
  expect_error(
    fit_partial_sqrt_lasso_shrinkage(
      X,
      y,
      core_set = 1:3,
      tested_set = 5:50,
      assume_independent = TRUE
    ),
    "partition"
  )
  expect_error(
    fit_partial_sqrt_lasso_shrinkage(
      X,
      y,
      core_set = 1:30,
      tested_set = 31:50,
      assume_independent = TRUE,
      bootstrap_B = 99L
    ),
    "strictly smaller"
  )
  expect_error(
    fit_partial_sqrt_lasso_shrinkage(
      X,
      y,
      core_set = 1:3,
      threshold_eta = 47,
      assume_independent = TRUE,
      bootstrap_B = 99L
    ),
    "strictly between zero and q"
  )
})

test_that("eta changes the endpoint but not the unthresholded max test", {
  set.seed(75)
  n <- 45L
  p <- 90L
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  multipliers <- matrix(rnorm(n * 99L), n, 99L)

  strict <- fit_partial_sqrt_lasso_shrinkage(
    X,
    y,
    core_set = 1:5,
    threshold_eta = 0.25,
    assume_independent = TRUE,
    bootstrap_multipliers = multipliers,
    max_iter = 10000L,
    tol = 1e-4
  )
  liberal <- fit_partial_sqrt_lasso_shrinkage(
    X,
    y,
    core_set = 1:5,
    threshold_eta = 10,
    assume_independent = TRUE,
    bootstrap_multipliers = multipliers,
    max_iter = 10000L,
    tol = 1e-4,
    null_geometry = strict$null_geometry,
    full_operator_norm = strict$full_solver$operator_norm
  )

  expect_equal(
    strict$inference$theta_tilde,
    liberal$inference$theta_tilde,
    tolerance = 1e-12
  )
  expect_equal(
    strict$inference$bootstrap_statistics,
    liberal$inference$bootstrap_statistics,
    tolerance = 1e-12
  )
  expect_gte(
    length(liberal$full_solver$support_tested),
    length(strict$full_solver$support_tested)
  )
})
