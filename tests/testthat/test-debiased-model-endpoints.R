test_that("submodel endpoints project the all-X full estimates when q exceeds n", {
  set.seed(40)
  n <- 32L
  p <- 70L
  p1 <- 6L
  tested <- (p1 + 1L):p
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(c(1.2, -1, 0.8, -0.7, 0.6, -0.5), rep(0, p - p1))
  y <- drop(X %*% beta + rnorm(n))

  fit <- fit_hd_shrinkage(
    X,
    y,
    restricted_set = tested,
    endpoint = "debiased_models",
    test = "max",
    precision = "diagonal",
    assume_independent = TRUE,
    bootstrap_B = 99L,
    bootstrap_seed = 41L,
    max_iter = 10000L,
    tol = 1e-4
  )

  expect_identical(fit$endpoint, "debiased_models")
  expect_gt(fit$dimensions[["q"]], fit$dimensions[["n"]])
  expect_equal(fit$beta$full, fit$beta$full_debiased)
  expect_equal(fit$beta$submodel, fit$beta$submodel_debiased)
  expect_equal(fit$beta$restricted, fit$beta$submodel)
  expect_equal(fit$beta$submodel[tested], rep(0, length(tested)))
  expect_equal(
    fit$beta$submodel[seq_len(p1)],
    fit$beta$full[seq_len(p1)],
    tolerance = 1e-12
  )
  expect_equal(
    fit$beta$submodel_regularized[seq_len(p1)],
    fit$beta$full_regularized[seq_len(p1)],
    tolerance = 1e-12
  )
  expect_equal(
    fit$beta$submodel_regularized[tested],
    rep(0, length(tested))
  )
  expect_equal(
    fit$beta$full[tested], fit$inference$theta_tilde,
    tolerance = 1e-12
  )
  expect_identical(fit$restriction$construction, "projection_of_full")
  expect_identical(fit$restricted_solver$construction, "projection_of_full")
  expect_identical(fit$restricted_solver$iterations, 0L)
  expect_true(fit$endpoint_precision$shared_full_fit)
  expect_identical(
    fit$endpoint_precision$submodel,
    fit$endpoint_precision$full
  )
  expect_gt(
    max(abs(fit$beta$full - fit$beta$full_regularized)),
    1e-8
  )
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
  expect_length(coef(fit, method = "FM"), p + 1L)
  expect_length(coef(fit, method = "SM"), p + 1L)
  expect_length(coef(fit, method = "FSL"), p + 1L)
  expect_length(coef(fit, method = "SMSL"), p + 1L)
  expect_length(coef(fit, method = "FDSL"), p + 1L)
  expect_length(coef(fit, method = "SMDSL"), p + 1L)
  expect_length(coef(fit, method = "PPS"), p + 1L)
})

test_that("projected debiasing rejects unsupported nonzero restrictions", {
  set.seed(42)
  X <- matrix(rnorm(25 * 40), 25, 40)
  y <- rnorm(25)
  expect_error(
    fit_hd_shrinkage(
      X,
      y,
      restricted_set = 31:40,
      t = rep(0.1, 10),
      endpoint = "debiased_models",
      precision = "diagonal",
      assume_independent = TRUE
    ),
    "coordinate restriction"
  )
})
