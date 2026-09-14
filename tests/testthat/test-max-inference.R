test_that("compiled multiplier maxima match the direct R calculation", {
  set.seed(20)
  n <- 12L
  q <- 17L
  B <- 103L
  psi <- matrix(rnorm(n * q), n, q)
  psi <- sweep(psi, 2L, colMeans(psi), "-")
  standard_error <- sqrt(colMeans(psi^2))
  multipliers <- matrix(rnorm(n * B), n, B)

  observed_cpp <- HDMaxShrink:::cpp_multiplier_max(
    psi,
    standard_error,
    multipliers,
    block_size = 19L
  )
  direct <- apply(
    crossprod(psi, multipliers) / sqrt(n) / standard_error,
    2L,
    function(value) max(abs(value))
  )
  expect_equal(as.numeric(observed_cpp), as.numeric(direct), tolerance = 1e-12)
})

test_that("maximum-test seed does not alter the caller RNG state", {
  set.seed(21)
  psi <- matrix(rnorm(20 * 8), 20, 8)
  psi <- sweep(psi, 2L, colMeans(psi), "-")
  before <- .Random.seed
  result <- max_test_hd(
    theta_tilde = rep(0, 8),
    target = rep(0, 8),
    psi_centered = psi,
    bootstrap_B = 99L,
    bootstrap_seed = 123L
  )
  expect_identical(.Random.seed, before)
  expect_equal(result$statistic, 0)
  expect_true(is.finite(result$null_energy))
  expect_gt(result$null_energy, 0)
})

test_that("max positive-part weight is bounded", {
  result <- HDMaxShrink:::.make_max_shrinkage(
    beta_full = c(2, -1),
    beta_restricted = c(0, 0),
    statistic = 4,
    null_energy = 5,
    reject = TRUE,
    epsilon = 1e-10
  )
  expect_gte(result$positive_weight, 0)
  expect_lte(result$positive_weight, 1)
  expect_equal(result$positive_part, result$positive_weight * c(2, -1))
  expect_equal(result$protected_weight, 1)
  expect_equal(result$protected_positive_part, c(2, -1))
})

test_that("protected positive-part equals PS when the max test does not reject", {
  result <- HDMaxShrink:::.make_max_shrinkage(
    beta_full = c(2, -1),
    beta_restricted = c(0.5, 0.25),
    statistic = 4,
    null_energy = 5,
    reject = FALSE,
    epsilon = 1e-10
  )
  expect_equal(result$protected_weight, result$positive_weight, tolerance = 0)
  expect_equal(
    result$protected_positive_part,
    result$positive_part,
    tolerance = 0
  )
  expect_gte(result$protected_weight, 0)
  expect_lte(result$protected_weight, 1)
})
