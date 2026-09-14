test_that("strict CPSS-LASSO is reproducible and exposes its audit trail", {
  set.seed(202608251)
  n <- 60L
  p <- 30L
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(2.2, -1.9, 1.6, rep(0, p - 3L))
  y <- drop(X %*% beta + rnorm(n, sd = 0.4))

  first <- cpss_select_core(
    X, y, selector = "lasso", complementary_pairs = 3L,
    base_selection_size = 5L, stability_threshold = 0.6,
    seed = 421L, path_points = 30L
  )
  second <- cpss_select_core(
    X, y, selector = "lasso", complementary_pairs = 3L,
    base_selection_size = 5L, stability_threshold = 0.6,
    seed = 421L, path_points = 30L
  )

  expect_s3_class(first, "cpss_core_selection")
  expect_identical(first$core_set, second$core_set)
  expect_equal(first$selection_frequency, second$selection_frequency,
               tolerance = 0)
  expect_length(first$selection_frequency, p)
  expect_true(all(first$selection_frequency >= 0 &
                    first$selection_frequency <= 1))
  expect_equal(nrow(first$base_fit_diagnostics), 6L)
  expect_true("half_sample_index" %in% names(first$base_fit_diagnostics))
  expect_identical(
    first$base_fit_diagnostics$half_sample_index,
    second$base_fit_diagnostics$half_sample_index
  )
  expect_true(all(first$base_fit_diagnostics$selected_count <= 5L))
  expect_true(first$mb_pfer_applicable)
  expect_match(first$pfer_bound_label, "requires_assumptions")
  expect_identical(first$path_lower_bound_semantics,
                   "glmnet_lambda_min_ratio")
  expect_true(all(1:3 %in% first$core_set))
  expect_true(first$no_top_k_fallback)
})

test_that("mandatory-core CPSS conditions on required variables", {
  set.seed(202608255)
  n <- 80L
  p <- 22L
  feature_names <- paste0("v", seq_len(p))
  mandatory_design <- matrix(rnorm(n * 2L), n, 2L)
  optional_base <- matrix(rnorm(n * (p - 3L)), n, p - 3L)
  aliased_optional <- mandatory_design[, 1L] -
    0.75 * mandatory_design[, 2L]
  X <- cbind(mandatory_design, aliased_optional, optional_base)
  colnames(X) <- feature_names
  y <- drop(
    1.8 * X[, 1L] - 1.4 * X[, 2L] +
      2.1 * X[, 4L] + rnorm(n, sd = 0.45)
  )

  selection <- cpss_select_core(
    X, y, selector = "lasso", complementary_pairs = 3L,
    base_selection_size = 4L, stability_threshold = 0.6,
    seed = 424L, path_points = 30L,
    mandatory_core = c("v1", "v2"),
    candidate_set = paste0("v", 3:20)
  )

  expect_true(selection$conditional_on_mandatory_core)
  expect_identical(selection$partialization,
                   "within_half_centered_FWL_economy_SVD")
  expect_identical(selection$mandatory_feature_names, c("v1", "v2"))
  expect_true(all(selection$mandatory_core %in% selection$core_set))
  expect_length(intersect(
    selection$mandatory_core, selection$selected_extension
  ), 0L)
  expect_true(3L %in% selection$ineligible_candidate_set)
  expect_true(all(is.na(
    selection$cpss_selection_frequency[selection$mandatory_core]
  )))
  expect_true(all(is.na(
    selection$selection_frequency[selection$mandatory_core]
  )))
  expect_equal(
    selection$retention_frequency[selection$mandatory_core],
    rep(1, 2L), tolerance = 0
  )
  expect_true(all(
    selection$base_fit_diagnostics$mandatory_rank == 2L
  ))
  expect_true(all(
    selection$base_fit_diagnostics$warning_count == 0L
  ))
  expect_identical(selection$pfer_family,
                   "eligible_optional_candidates_only")
  expect_equal(
    selection$pfer_upper_bound_mb,
    4^2 / ((2 * 0.6 - 1) *
      length(selection$eligible_candidate_set))
  )
  expect_true(all(
    selection$stability_table$core_role[
      selection$stability_table$feature %in% c("v1", "v2")
    ] == "mandatory"
  ))
  expect_true(all(
    selection$stability_table$core_role[
      selection$stability_table$feature %in% c("v21", "v22")
    ] == "outside_candidate_universe"
  ))
})

test_that("mandatory-core CPSS is invariant to mandatory linear components", {
  set.seed(202608256)
  n <- 76L
  p0 <- 3L
  q <- 16L
  X0 <- matrix(rnorm(n * p0), n, p0)
  C <- matrix(rnorm(n * q), n, q)
  loading <- matrix(rnorm(p0 * q), p0, q)
  X_plain <- cbind(X0, C)
  X_shifted <- cbind(X0, C + X0 %*% loading)
  feature_names <- c(paste0("core", seq_len(p0)),
                     paste0("gene", seq_len(q)))
  colnames(X_plain) <- colnames(X_shifted) <- feature_names
  y <- drop(X0 %*% c(1.2, -0.9, 0.7) +
              C[, 1L] * 1.8 - C[, 2L] * 1.4 + rnorm(n, sd = 0.5))
  arguments <- list(
    y = y, selector = "lasso", complementary_pairs = 3L,
    base_selection_size = 4L, stability_threshold = 0.6,
    seed = 425L, path_points = 35L,
    mandatory_core = paste0("core", seq_len(p0))
  )

  plain <- do.call(cpss_select_core, c(list(X = X_plain), arguments))
  shifted <- do.call(cpss_select_core, c(list(X = X_shifted), arguments))

  expect_equal(
    plain$cpss_selection_frequency,
    shifted$cpss_selection_frequency,
    tolerance = 1e-10
  )
  expect_identical(
    plain$selected_extension_feature_names,
    shifted$selected_extension_feature_names
  )
})

test_that("CPSS feature identities align a reordered analysis design", {
  set.seed(202608252)
  p <- 24L
  feature_names <- paste0("gene_", seq_len(p))
  X_selection <- matrix(rnorm(70L * p), 70L, p,
                        dimnames = list(NULL, feature_names))
  beta <- c(2.4, -2, 1.8, rep(0, p - 3L))
  y_selection <- drop(X_selection %*% beta + rnorm(70L, sd = 0.35))
  selection <- cpss_select_core(
    X_selection, y_selection, selector = "lasso",
    complementary_pairs = 3L, base_selection_size = 5L,
    stability_threshold = 0.6, seed = 422L, path_points = 30L
  )

  X_analysis <- matrix(rnorm(45L * p), 45L, p,
                       dimnames = list(NULL, feature_names))
  y_analysis <- drop(X_analysis %*% beta + rnorm(45L, sd = 0.6))
  permutation <- sample(seq_len(p))
  fit <- fit_cpss_sm(
    X_analysis[, permutation, drop = FALSE], y_analysis,
    selection = selection
  )

  expect_s3_class(fit, "cpss_sm_fit")
  expect_setequal(
    fit$feature_names[fit$core_set], selection$core_feature_names
  )
  expect_equal(unname(fit$beta[fit$tested_set]),
               numeric(length(fit$tested_set)),
               tolerance = 0)
  expect_equal(fit$solver$restriction_violation, 0, tolerance = 0)
  expect_equal(fit$solver$rank, length(fit$core_set))
  prediction <- predict(
    fit, X_analysis[1:4, rev(seq_len(p)), drop = FALSE]
  )
  expect_length(prediction, 4L)
  expect_length(coef(fit), p + 1L)

  gaussian_draws <- matrix(rnorm(45L * 99L), 45L, 99L)
  shrinkage_fit <- fit_cpss_ridge_shrinkage(
    X_analysis[, permutation, drop = FALSE], y_analysis,
    selection = selection, ridge_lambda = 0.4,
    gaussian_draws = gaussian_draws, bootstrap_B = 99L
  )
  expect_setequal(
    shrinkage_fit$feature_names[shrinkage_fit$core_set],
    selection$core_feature_names
  )
  expect_identical(
    shrinkage_fit$inference$selection_conditioning,
    "conditional on an independently learned CPSS core"
  )
})

test_that("CPSS--SM/Ridge uses all X and preserves estimator identities", {
  set.seed(202608253)
  n <- 46L
  p <- 72L
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(1.4, -1.2, 1, 0.8, rep(0, p - 4L))
  y <- drop(X %*% beta + rnorm(n, sd = 0.8))
  gaussian_draws <- matrix(rnorm(n * 99L), n, 99L)
  lambda <- 0.35

  fit <- fit_cpss_ridge_shrinkage(
    X, y, core_set = 1:4, ridge_lambda = lambda,
    selector_label = "fixed-test-core", gaussian_draws = gaussian_draws,
    bootstrap_B = 99L, shrinkage_calibration = "inverse_moment"
  )

  expect_s3_class(fit, "hd_shrinkage_fit")
  expect_identical(fit$endpoint, "cpss_ridge_exact_null")
  expect_identical(fit$full_solver$construction, "all_X_dual_Ridge")
  expect_identical(fit$test, "max_partial_t")
  expect_identical(fit$test_calibration, "conditional_gaussian")
  expect_identical(
    fit$inference$shrinkage_calibration_type, "inverse_moment"
  )
  expect_equal(fit$beta$submodel[5:p], numeric(p - 4L), tolerance = 0)
  expect_equal(fit$restricted_solver$restriction_violation, 0, tolerance = 0)
  expect_equal(fit$inference$q_original, p - 4L)
  expect_lte(fit$inference$q_effective, p - 4L)
  expect_gt(fit$inference$q_effective, 0L)
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

  prep <- HDMaxShrink:::.prepare_design(X, y, scale_y = TRUE)
  direct_standardized <- drop(t(prep$X) %*% solve(
    prep$X %*% t(prep$X) + n * lambda * diag(n), prep$y
  ))
  direct <- direct_standardized * prep$y_scale / prep$x_scale
  expect_equal(fit$beta$full, direct, tolerance = 1e-10)
  expect_lt(fit$full_solver$max_kkt_error, 1e-9)
  expect_true(all(unlist(fit$standardized_intercept) == 0))
  expect_equal(unname(coef(fit, method = "FM")[-1L]), fit$beta$full)
  expect_equal(unname(coef(fit, method = "SM")[-1L]), fit$beta$submodel)
})

test_that("the optional CPSS-MCP selector records grpreg path diagnostics", {
  skip_if_not_installed("grpreg")
  set.seed(202608254)
  X <- matrix(rnorm(54L * 25L), 54L, 25L)
  y <- drop(X[, 1L] * 2 - X[, 2L] * 1.5 + rnorm(54L, sd = 0.5))
  selection <- cpss_select_core(
    X, y, selector = "mcp", complementary_pairs = 1L,
    base_selection_size = 4L, stability_threshold = 0.5,
    seed = 423L, path_points = 20L, mcp_max_iter = 10000L
  )
  repeated <- cpss_select_core(
    X, y, selector = "mcp", complementary_pairs = 1L,
    base_selection_size = 4L, stability_threshold = 0.5,
    seed = 423L, path_points = 20L, mcp_max_iter = 10000L
  )
  expect_identical(selection$base_selector, "mcp")
  expect_identical(selection$core_set, repeated$core_set)
  expect_equal(selection$selection_frequency,
               repeated$selection_frequency, tolerance = 0)
  expect_true(
    "selected_path_locally_convex" %in%
      names(selection$base_fit_diagnostics)
  )
  expect_false(selection$mb_pfer_applicable)
  expect_true(is.na(selection$pfer_upper_bound_mb))
  expect_identical(selection$path_lower_bound_semantics,
                   "grpreg_lambda_min_ratio")
  expect_identical(selection$base_selector_engine,
                   "grpreg_singleton_group_MCP")
  expect_identical(selection$mcp_local_convexity_diagnostic,
                   "not_available_from_grpreg")
  expect_true(all(is.na(
    selection$base_fit_diagnostics$selected_path_locally_convex
  )))
  expect_true(all(selection$base_fit_diagnostics$converged))
  expect_true(all(selection$base_fit_diagnostics$warning_count == 0L))
  expect_lte(max(selection$base_fit_diagnostics$selected_count), 4L)
})
