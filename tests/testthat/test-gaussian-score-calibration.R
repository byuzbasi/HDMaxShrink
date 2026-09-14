test_that("conditional Gaussian score kernel matches direct formulas", {
  set.seed(81)
  n <- 32L
  p <- 48L
  core <- 1:4
  tested <- 5:p
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  prep <- HDMaxShrink:::.prepare_design(X, y)
  geometry <- HDMaxShrink:::.make_partial_null_geometry(
    X = prep$X,
    core_set = core,
    tested_set = tested,
    rank_tolerance = 1e-10
  )
  null_state <- HDMaxShrink:::cpp_partial_null_apply(
    prep$y,
    geometry$coefficient_operator,
    geometry$orthonormal_basis,
    geometry$residualized_tested
  )
  gaussian_draws <- matrix(rnorm(n * 199L), n, 199L)

  observed <- HDMaxShrink:::.gaussian_score_max_inference(
    null_state = null_state,
    geometry = geometry,
    alpha = 0.05,
    gaussian_draws = gaussian_draws,
    bootstrap_block_size = 31L
  )

  V <- geometry$residualized_tested
  Q <- geometry$orthonormal_basis
  residual <- null_state$residual
  residual_df <- n - ncol(Q) - 1L
  direction_norm <- sqrt(colSums(V^2))
  sigma_hat <- sqrt(sum(residual^2) / residual_df)
  expected_score <- drop(crossprod(V, residual)) /
    direction_norm / sigma_hat

  residual_draws <- sweep(gaussian_draws, 2L, colMeans(gaussian_draws), "-")
  residual_draws <- residual_draws - Q %*% crossprod(Q, residual_draws)
  draw_scale <- sqrt(colSums(residual_draws^2) / residual_df)
  standardized <- sweep(crossprod(V, residual_draws), 1L, direction_norm, "/")
  standardized <- sweep(standardized, 2L, draw_scale, "/")
  expected_bootstrap <- apply(abs(standardized), 2L, max)

  expect_equal(observed$residual_df, as.numeric(residual_df))
  expect_equal(observed$sigma_hat, sigma_hat, tolerance = 1e-12)
  expect_equal(observed$direction_norm, direction_norm, tolerance = 1e-12)
  expect_equal(observed$score, expected_score, tolerance = 1e-11)
  expect_equal(observed$statistic, max(abs(expected_score)), tolerance = 1e-11)
  expect_equal(
    observed$bootstrap_statistics,
    expected_bootstrap,
    tolerance = 1e-11
  )
})

test_that("Monte Carlo decision uses an exact finite-draw p-value", {
  bootstrap_statistics <- as.numeric(1:99)
  tied <- HDMaxShrink:::.monte_carlo_max_decision(
    observed = 95,
    bootstrap_statistics = bootstrap_statistics,
    alpha = 0.05
  )
  above <- HDMaxShrink:::.monte_carlo_max_decision(
    observed = 95.1,
    bootstrap_statistics = bootstrap_statistics,
    alpha = 0.05
  )

  expect_equal(tied$critical_value, 95)
  expect_equal(tied$p_value, 0.06)
  expect_false(tied$reject)
  expect_equal(above$p_value, 0.05)
  expect_true(above$reject)
})

test_that("conditional Gaussian score calibration is scale invariant", {
  set.seed(82)
  n <- 28L
  p <- 50L
  core <- 1:3
  tested <- 4:p
  X <- matrix(rnorm(n * p), n, p)
  prep <- HDMaxShrink:::.prepare_design(X, rnorm(n))
  geometry <- HDMaxShrink:::.make_partial_null_geometry(
    X = prep$X,
    core_set = core,
    tested_set = tested,
    rank_tolerance = 1e-10
  )
  null_state <- HDMaxShrink:::cpp_partial_null_apply(
    prep$y,
    geometry$coefficient_operator,
    geometry$orthonormal_basis,
    geometry$residualized_tested
  )
  gaussian_draws <- matrix(rnorm(n * 99L), n, 99L)
  baseline <- HDMaxShrink:::.gaussian_score_max_inference(
    null_state,
    geometry,
    gaussian_draws = gaussian_draws
  )
  rescaled_state <- null_state
  rescaled_state$residual <- 3.7 * rescaled_state$residual
  rescaled <- HDMaxShrink:::.gaussian_score_max_inference(
    rescaled_state,
    geometry,
    gaussian_draws = gaussian_draws
  )

  expect_equal(rescaled$statistic, baseline$statistic, tolerance = 1e-12)
  expect_equal(
    rescaled$bootstrap_statistics,
    baseline$bootstrap_statistics,
    tolerance = 0
  )
  expect_equal(rescaled$p_value, baseline$p_value)
  expect_identical(rescaled$reject, baseline$reject)
})

test_that("score maxima transform monotonically to max partial-t", {
  set.seed(83)
  n <- 36L
  p <- 65L
  core <- 1:5
  tested <- 6:p
  X <- matrix(rnorm(n * p), n, p)
  prep <- HDMaxShrink:::.prepare_design(X, rnorm(n))
  geometry <- HDMaxShrink:::.make_partial_null_geometry(
    X = prep$X,
    core_set = core,
    tested_set = tested,
    rank_tolerance = 1e-10
  )
  null_state <- HDMaxShrink:::cpp_partial_null_apply(
    prep$y,
    geometry$coefficient_operator,
    geometry$orthonormal_basis,
    geometry$residualized_tested
  )
  gaussian_draws <- matrix(rnorm(n * 199L), n, 199L)
  score <- HDMaxShrink:::.gaussian_score_max_inference(
    null_state,
    geometry,
    gaussian_draws = gaussian_draws
  )
  partial_t <- HDMaxShrink:::.as_max_partial_t_inference(score)
  partial_t <- HDMaxShrink:::.set_max_shrinkage_calibration(
    partial_t, "inverse_moment"
  )

  expected <- sqrt(
    (score$residual_df - 1) * score$statistic^2 /
      (score$residual_df - score$statistic^2)
  )
  expected_inverse <- 1 / mean(1 / partial_t$bootstrap_statistics^2)

  expect_identical(partial_t$type, "gaussian_max_partial_t")
  expect_identical(
    partial_t$subtype,
    "conditional_gaussian_exact_null_max_partial_t"
  )
  expect_equal(partial_t$score_statistic, score$statistic, tolerance = 0)
  expect_equal(partial_t$max_partial_t, expected, tolerance = 1e-13)
  expect_equal(partial_t$statistic, partial_t$max_partial_t, tolerance = 0)
  expect_equal(partial_t$p_value, score$p_value, tolerance = 0)
  expect_identical(partial_t$reject, score$reject)
  expect_equal(
    partial_t$null_second_moment,
    mean(partial_t$bootstrap_statistics^2),
    tolerance = 1e-14
  )
  expect_equal(
    partial_t$null_inverse_moment_calibration,
    expected_inverse,
    tolerance = 1e-14
  )
  expect_equal(partial_t$shrinkage_calibration, expected_inverse)
  expect_identical(partial_t$shrinkage_calibration_type, "inverse_moment")
  expect_equal(partial_t$null_energy, partial_t$null_second_moment)
})

test_that("conditional Gaussian calibration requires at least 99 draws", {
  set.seed(84)
  n <- 25L
  X <- matrix(rnorm(n * 40L), n, 40L)
  prep <- HDMaxShrink:::.prepare_design(X, rnorm(n))
  geometry <- HDMaxShrink:::.make_partial_null_geometry(
    prep$X, 1:3, 4:40, 1e-10
  )
  null_state <- HDMaxShrink:::cpp_partial_null_apply(
    prep$y,
    geometry$coefficient_operator,
    geometry$orthonormal_basis,
    geometry$residualized_tested
  )
  expect_error(
    HDMaxShrink:::.gaussian_score_max_inference(
      null_state,
      geometry,
      gaussian_draws = matrix(rnorm(n * 98L), n, 98L)
    ),
    "at least 99"
  )
})
