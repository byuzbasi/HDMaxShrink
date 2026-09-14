.prepare_cpss_null_endpoint <- function(
    prep, core_set, rank_tolerance = 1e-10) {
  n <- nrow(prep$X)
  p <- ncol(prep$X)
  tested_set <- setdiff(seq_len(p), core_set)
  if (!length(tested_set)) {
    stop("The CPSS--SM core must leave at least one tested coordinate.",
         call. = FALSE)
  }
  rank_tolerance <- .check_positive_scalar(
    rank_tolerance, "rank_tolerance"
  )

  core_design <- prep$X[, core_set, drop = FALSE]
  core_svd <- svd(
    core_design, nu = length(core_set), nv = length(core_set)
  )
  cutoff <- rank_tolerance * max(dim(core_design)) * core_svd$d[1L]
  rank <- sum(core_svd$d > cutoff)
  if (rank != length(core_set)) {
    stop(
      "The selected core is not full column rank at rank_tolerance.",
      call. = FALSE
    )
  }
  core_basis <- core_svd$u[, seq_len(rank), drop = FALSE]
  tested_design <- prep$X[, tested_set, drop = FALSE]
  residualized <- tested_design -
    core_basis %*% crossprod(core_basis, tested_design)
  original_norm <- sqrt(colSums(tested_design^2))
  residualized_norm <- sqrt(colSums(residualized^2))
  relative_norm <- residualized_norm / pmax(
    original_norm, sqrt(.Machine$double.eps)
  )
  retained <- is.finite(relative_norm) &
    relative_norm > sqrt(.Machine$double.eps)
  if (!any(retained)) {
    stop(
      "No testable coordinate remains after residualizing against the core.",
      call. = FALSE
    )
  }
  effective_tested_set <- tested_set[retained]
  geometry <- .make_partial_null_geometry(
    X = prep$X,
    core_set = core_set,
    tested_set = effective_tested_set,
    rank_tolerance = rank_tolerance
  )
  null_state <- cpp_partial_null_apply(
    y = prep$y,
    coefficient_operator = geometry$coefficient_operator,
    orthonormal_basis = geometry$orthonormal_basis,
    residualized_tested = geometry$residualized_tested
  )
  beta_standardized <- numeric(p)
  beta_standardized[core_set] <- null_state$beta_core
  coefficient_backscale <- prep$y_scale / prep$x_scale
  beta <- beta_standardized * coefficient_backscale
  intercept <- prep$y_center - sum(prep$x_center * beta)
  dropped_set <- tested_set[!retained]
  list(
    beta = beta,
    intercept = intercept,
    beta_standardized = beta_standardized,
    core_set = core_set,
    tested_set = tested_set,
    effective_tested_set = effective_tested_set,
    dropped_tested_set = dropped_set,
    q_original = length(tested_set),
    q_effective = length(effective_tested_set),
    geometry = geometry,
    null_state = null_state,
    rank = geometry$rank,
    condition_number = geometry$condition_number,
    relative_tested_norm = relative_norm,
    rank_tolerance = rank_tolerance
  )
}

#' Fit the CPSS-guided exact-null submodel endpoint
#'
#' The selected core is refitted from scratch using an economy-SVD
#' Moore--Penrose operator. All coefficients outside the core are set exactly
#' to zero. The function never copies full-model coefficients and never forms
#' an OLS Gram inverse.
#'
#' @param X Numeric analysis-sample design matrix.
#' @param y Numeric analysis-sample response.
#' @param core_set One-based selected core indices, or a
#'   `cpss_core_selection` object.
#' @param selection Alternatively, an object returned by
#'   [cpss_select_core()]. Named features are aligned to `X` by identity.
#' @param rank_tolerance Relative economy-SVD rank tolerance.
#'
#' @return A `cpss_sm_fit` object with original-scale coefficients and
#'   exact-null diagnostics.
#' @references
#' Belloni, A. and Chernozhukov, V. (2013). Least squares after model
#' selection in high-dimensional sparse models. *Bernoulli*, 19, 521--547.
#' \doi{10.3150/11-BEJ410}
#' @export
fit_cpss_sm <- function(
    X,
    y,
    core_set = NULL,
    selection = NULL,
    rank_tolerance = 1e-10) {
  call <- match.call()
  X <- .validate_numeric_matrix(X, "X")
  y <- .validate_response(y, nrow(X))
  resolved <- .resolve_cpss_core(X, core_set, selection)
  prep <- .prepare_design(X, y, scale_y = TRUE)
  state <- .prepare_cpss_null_endpoint(
    prep, resolved$core_set, rank_tolerance
  )
  feature_state <- .cpss_feature_names(X)
  names(state$beta) <- feature_state$names
  output <- list(
    call = call,
    coefficients = state$beta,
    beta = state$beta,
    intercept = as.numeric(state$intercept),
    standardized_intercept = 0,
    core_set = state$core_set,
    tested_set = state$tested_set,
    effective_tested_set = state$effective_tested_set,
    dropped_tested_set = state$dropped_tested_set,
    dimensions = c(
      n = nrow(X), p = ncol(X), p1 = length(state$core_set),
      q = state$q_original, q_effective = state$q_effective
    ),
    solver = list(
      construction = "CPSS_guided_exact_null_refit_svd",
      rank = state$rank,
      condition_number = state$condition_number,
      rank_tolerance = state$rank_tolerance,
      residual_norm = state$null_state$residual_norm,
      converged = TRUE,
      restriction_violation = max(abs(state$beta[state$tested_set]))
    ),
    restriction = list(
      type = "coordinate_exact_null",
      core_set = state$core_set,
      tested_set = state$tested_set,
      target = numeric(state$q_original)
    ),
    selection = resolved$selection,
    feature_names = feature_state$names,
    preprocessing = prep[c(
      "x_center", "x_scale", "y_center", "y_scale",
      "standardized_intercept"
    )]
  )
  class(output) <- c("cpss_sm_fit", "list")
  output
}

#' @export
coef.cpss_sm_fit <- function(object, ...) {
  result <- c(`(Intercept)` = object$intercept, object$beta)
  names(result)[-1L] <- object$feature_names
  result
}

#' @export
predict.cpss_sm_fit <- function(object, newdata, ...) {
  newdata <- .validate_numeric_matrix(newdata, "newdata")
  if (ncol(newdata) != length(object$beta)) {
    stop("newdata must have the original feature universe.", call. = FALSE)
  }
  if (!is.null(colnames(newdata))) {
    current <- .cpss_feature_names(newdata)$names
    if (!setequal(current, object$feature_names)) {
      stop("newdata feature names do not match the fitted model.",
           call. = FALSE)
    }
    newdata <- newdata[, match(object$feature_names, current), drop = FALSE]
  }
  as.numeric(object$intercept + newdata %*% object$beta)
}

#' CPSS--SM shrinkage toward an all-predictor dual-Ridge endpoint
#'
#' `fit_cpss_ridge_shrinkage()` combines an all-`X` Ridge full-model endpoint
#' with a fresh exact-null SVD submodel endpoint. A conditional-Gaussian
#' maximum partial-t test is calculated from the exact-null residual, not
#' from the biased Ridge estimator. The function returns preliminary-test,
#' Stein-type, positive-part, and pretest-protected positive-part estimators.
#'
#' @inheritParams fit_cpss_sm
#' @param ridge_lambda Optional fixed positive Ridge penalty. If omitted, GCV
#'   selects from `ridge_lambda_grid` using the existing all-`X` dual path.
#' @param ridge_lambda_grid Optional positive GCV grid.
#' @param selector_label Label describing the independently selected core.
#' @param alpha Conditional Monte Carlo test level.
#' @param gaussian_draws Optional reusable `n` by `bootstrap_B` matrix of iid
#'   standard Gaussian draws.
#' @param bootstrap_B Number of Gaussian draws when `gaussian_draws` is absent.
#' @param bootstrap_seed Optional reproducible Gaussian-draw seed.
#' @param bootstrap_block_size C++ calibration block size.
#' @param shrinkage_calibration Null calibration for the max-Stein weight.
#' @param epsilon Positive numerical floor for the squared statistic.
#'
#' @return An `hd_shrinkage_fit` object with `FM`, `SM`, `PT`, `S`, `PS`, and
#'   `PPS` endpoints.
#' @references
#' Shah, R. D. and Samworth, R. J. (2013). Variable selection with error
#' control: another look at stability selection. *Journal of the Royal
#' Statistical Society: Series B*, 75, 55--80.
#' \doi{10.1111/j.1467-9868.2011.01034.x}
#'
#' Dufour, J.-M. (2006). Monte Carlo tests with nuisance parameters: a general
#' approach to finite-sample inference and nonstandard asymptotics.
#' *Journal of Econometrics*, 133, 443--477.
#' \doi{10.1016/j.jeconom.2005.06.007}
#'
#' Yuzbasi, B., Arashi, M. and Ahmed, S. E. (2020). Shrinkage estimation
#' strategies in high-dimensional linear models. *International Statistical
#' Review*. \doi{10.1111/insr.12351}
#' @export
fit_cpss_ridge_shrinkage <- function(
    X,
    y,
    core_set = NULL,
    selection = NULL,
    ridge_lambda = NULL,
    ridge_lambda_grid = NULL,
    selector_label = NULL,
    alpha = 0.05,
    gaussian_draws = NULL,
    bootstrap_B = 999L,
    bootstrap_seed = NULL,
    bootstrap_block_size = 256L,
    rank_tolerance = 1e-10,
    shrinkage_calibration = c("inverse_moment", "second_moment"),
    epsilon = 1e-10) {
  call <- match.call()
  X <- .validate_numeric_matrix(X, "X")
  y <- .validate_response(y, nrow(X))
  resolved <- .resolve_cpss_core(X, core_set, selection)
  shrinkage_calibration <- match.arg(shrinkage_calibration)
  epsilon <- .check_positive_scalar(epsilon, "epsilon")
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie strictly between zero and one.", call. = FALSE)
  }
  if (is.null(selector_label)) {
    selector_label <- resolved$selection$selector %||% "fixed-core"
  }
  if (length(selector_label) != 1L || !is.character(selector_label) ||
      is.na(selector_label) || !nzchar(selector_label)) {
    stop("selector_label must be one nonempty string.", call. = FALSE)
  }

  prep <- .prepare_design(X, y, scale_y = TRUE)
  n <- nrow(prep$X)
  p <- ncol(prep$X)
  null <- .prepare_cpss_null_endpoint(
    prep, resolved$core_set, rank_tolerance
  )
  ridge <- .select_ridge_endpoint(
    prep$X, prep$y,
    lambda = ridge_lambda,
    lambda_grid = ridge_lambda_grid
  )
  beta_ridge_standardized <- as.numeric(ridge$beta)
  ridge_residual <- as.numeric(
    prep$X %*% beta_ridge_standardized - prep$y
  )
  ridge_gradient <- as.numeric(
    crossprod(prep$X, ridge_residual) / n +
      ridge$lambda * beta_ridge_standardized
  )
  beta_ridge <- beta_ridge_standardized *
    prep$y_scale / prep$x_scale
  beta_sm <- null$beta

  inference <- .gaussian_score_max_inference(
    null_state = null$null_state,
    geometry = null$geometry,
    alpha = alpha,
    bootstrap_B = bootstrap_B,
    bootstrap_seed = bootstrap_seed,
    gaussian_draws = gaussian_draws,
    bootstrap_block_size = bootstrap_block_size
  )
  inference <- .as_max_partial_t_inference(inference)
  inference <- .set_max_shrinkage_calibration(
    inference, shrinkage_calibration
  )
  inference$test_calibration <- "conditional_gaussian"
  inference$q_original <- null$q_original
  inference$q_effective <- null$q_effective
  inference$tested_set <- null$tested_set
  inference$effective_tested_set <- null$effective_tested_set
  inference$dropped_tested_set <- null$dropped_tested_set
  inference$selection_conditioning <- if (is.null(resolved$selection)) {
    "core treated as fixed before analysis"
  } else {
    "conditional on an independently learned CPSS core"
  }

  shrinkage <- .make_max_shrinkage(
    beta_full = beta_ridge,
    beta_restricted = beta_sm,
    statistic = inference$statistic,
    reject = inference$reject,
    epsilon = epsilon,
    shrinkage_calibration = inference$shrinkage_calibration,
    shrinkage_calibration_type =
      inference$shrinkage_calibration_type
  )
  beta <- list(
    full = beta_ridge,
    submodel = beta_sm,
    restricted = beta_sm,
    ridge = beta_ridge,
    mcp = NULL,
    full_regularized = beta_ridge,
    submodel_regularized = beta_sm,
    full_debiased = NULL,
    submodel_debiased = NULL,
    preliminary_test = shrinkage$preliminary_test,
    stein = shrinkage$stein,
    positive_part = shrinkage$positive_part,
    protected_positive_part = shrinkage$protected_positive_part
  )
  intercept <- lapply(beta, function(coefficient) {
    if (is.null(coefficient)) return(NULL)
    prep$y_center - sum(prep$x_center * coefficient)
  })
  standardized_intercept <- lapply(beta, function(coefficient) {
    if (is.null(coefficient)) return(NULL)
    0
  })
  feature_state <- .cpss_feature_names(X)
  ridge$beta <- NULL
  ridge$construction <- "all_X_dual_Ridge"
  ridge$penalty_scope <- "all_X"
  ridge$standardized_response <- TRUE
  ridge$effective_df <- ridge$degrees_freedom
  ridge$max_kkt_error <- max(abs(ridge_gradient))
  ridge$kkt_residual <- ridge$max_kkt_error
  restricted_solver <- list(
    beta = null$beta_standardized,
    beta_core = null$beta_standardized[null$core_set],
    beta_tested = null$beta_standardized[null$tested_set],
    converged = TRUE,
    iterations = 1L,
    objective = null$null_state$residual_norm / sqrt(n),
    restriction_violation = max(abs(beta_sm[null$tested_set])),
    construction = "CPSS_guided_exact_null_refit_svd",
    rank = null$rank,
    condition_number = null$condition_number,
    residual_norm = null$null_state$residual_norm,
    q_original = null$q_original,
    q_effective = null$q_effective,
    dropped_tested_set = null$dropped_tested_set
  )
  output <- list(
    call = call,
    dimensions = c(
      n = n, p = p, p1 = length(null$core_set),
      p2 = null$q_original, m = null$q_original, q = null$q_original
    ),
    core_set = null$core_set,
    tested_set = null$tested_set,
    selector_label = selector_label,
    selection = resolved$selection,
    feature_names = feature_state$names,
    beta = beta,
    intercept = intercept,
    standardized_intercept = standardized_intercept,
    endpoint = "cpss_ridge_exact_null",
    full_endpoint = "all_x_dual_ridge",
    test = "max_partial_t",
    test_calibration = "conditional_gaussian",
    reject = inference$reject,
    alpha = alpha,
    inference = inference,
    shrinkage = shrinkage,
    full_solver = ridge,
    restricted_solver = restricted_solver,
    restriction_projection = restricted_solver,
    restriction = list(
      type = "coordinate_exact_null",
      core_set = null$core_set,
      tested_set = null$tested_set,
      effective_tested_set = null$effective_tested_set,
      dropped_tested_set = null$dropped_tested_set,
      t = numeric(null$q_original),
      zero_block = TRUE,
      construction = "CPSS_guided_exact_null_refit_svd",
      exact_endpoint_restriction = TRUE,
      selection_conditioning = inference$selection_conditioning
    ),
    precision = list(method = "not_required", required = FALSE),
    score_representation =
      "conditional_gaussian_exact_null_partial_t",
    preprocessing = prep[c(
      "x_center", "x_scale", "y_center", "y_scale",
      "standardized_intercept"
    )],
    null_geometry = null$geometry,
    geometry_diagnostics = list(
      rank = null$rank,
      condition_number = null$condition_number,
      rank_tolerance = rank_tolerance,
      q_original = null$q_original,
      q_effective = null$q_effective,
      dropped_tested_set = null$dropped_tested_set
    )
  )
  class(output) <- "hd_shrinkage_fit"
  output
}
