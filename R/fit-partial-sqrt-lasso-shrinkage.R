.partial_design_signature <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  row_weight <- sin(seq_len(n)) + cos(seq_len(n) / 3)
  column_weight <- sin(seq_len(p) / 5) - cos(seq_len(p) / 7)
  c(
    n = n,
    p = p,
    sum_squares = sum(X^2),
    weighted_rows = sum(row_weight * rowSums(X)),
    weighted_columns = sum(column_weight * colSums(X^2))
  )
}

.partial_partition <- function(core_set, tested_set, p) {
  core_set <- .validate_indices(core_set, p, "core_set")
  tested_set <- tested_set %||% setdiff(seq_len(p), core_set)
  tested_set <- .validate_indices(tested_set, p, "tested_set")
  if (length(intersect(core_set, tested_set))) {
    stop("core_set and tested_set must be disjoint.", call. = FALSE)
  }
  if (!setequal(c(core_set, tested_set), seq_len(p))) {
    stop(
      "core_set and tested_set must form a partition of all p columns.",
      call. = FALSE
    )
  }
  list(core_set = core_set, tested_set = tested_set)
}

.make_partial_null_geometry <- function(
    X,
    core_set,
    tested_set,
    rank_tolerance) {
  geometry <- cpp_partial_null_geometry(
    X_core = X[, core_set, drop = FALSE],
    X_tested = X[, tested_set, drop = FALSE],
    rank_tolerance = rank_tolerance
  )
  geometry$n <- nrow(X)
  geometry$p <- ncol(X)
  geometry$p1 <- length(core_set)
  geometry$q <- length(tested_set)
  geometry$core_set <- as.integer(core_set)
  geometry$tested_set <- as.integer(tested_set)
  geometry$design_signature <- .partial_design_signature(X)
  class(geometry) <- c("hd_partial_null_geometry", "list")
  geometry
}

.validate_partial_null_geometry <- function(
    geometry,
    X,
    core_set,
    tested_set) {
  required <- c(
    "coefficient_operator", "orthonormal_basis", "residualized_tested",
    "residualized_second_moment", "rank", "condition_number",
    "orthogonality_error", "n", "p", "p1", "q", "core_set",
    "tested_set", "design_signature"
  )
  if (!inherits(geometry, "hd_partial_null_geometry") ||
      !all(required %in% names(geometry))) {
    stop(
      "null_geometry must be created by this partial-null implementation.",
      call. = FALSE
    )
  }
  expected_dimensions <- c(
    n = nrow(X), p = ncol(X), p1 = length(core_set), q = length(tested_set)
  )
  observed_dimensions <- unlist(geometry[names(expected_dimensions)])
  if (!isTRUE(all.equal(
    as.numeric(observed_dimensions), as.numeric(expected_dimensions),
    tolerance = 0
  )) ||
      !identical(as.integer(geometry$core_set), as.integer(core_set)) ||
      !identical(as.integer(geometry$tested_set), as.integer(tested_set))) {
    stop("null_geometry does not match the requested partition.", call. = FALSE)
  }
  signature <- .partial_design_signature(X)
  if (!isTRUE(all.equal(
    as.numeric(geometry$design_signature), as.numeric(signature),
    tolerance = 1e-12, scale = 1
  ))) {
    stop("null_geometry was computed from a different design.", call. = FALSE)
  }
  if (!is.matrix(geometry$residualized_tested) ||
      any(dim(geometry$residualized_tested) != c(nrow(X), length(tested_set))) ||
      any(!is.finite(geometry$residualized_tested))) {
    stop("null_geometry contains an invalid residualized design.", call. = FALSE)
  }
  invisible(geometry)
}

.partial_debiased_endpoint_state <- function(
    null_state,
    geometry,
    beta_tested_standardized,
    tested_scale) {
  score_state_standardized <- cpp_partial_debiased_scores(
    y_residualized = null_state$residual,
    residualized_tested = geometry$residualized_tested,
    beta_tested = beta_tested_standardized
  )
  theta_original <- as.numeric(
    score_state_standardized$theta_tilde
  ) / tested_scale
  psi_original <- sweep(
    score_state_standardized$psi_centered,
    2L,
    tested_scale,
    "/"
  )
  variance <- colMeans(psi_original^2)
  if (any(!is.finite(variance)) || any(variance <= 0)) {
    stop("Every partial-debiased score variance must be positive.", call. = FALSE)
  }
  list(
    theta_tilde = theta_original,
    standard_error = sqrt(variance),
    variance = variance,
    theta_tilde_standardized = as.numeric(
      score_state_standardized$theta_tilde
    ),
    correction_standardized = as.numeric(
      score_state_standardized$correction
    ),
    pilot_residual_norm = score_state_standardized$residual_norm,
    score_second_moment = score_state_standardized$second_moment,
    psi_centered = psi_original
  )
}

.attach_partial_endpoint_state <- function(inference, endpoint_state) {
  inference$theta_tilde <- endpoint_state$theta_tilde
  inference$standard_error <- endpoint_state$standard_error
  inference$variance <- endpoint_state$variance
  inference$theta_tilde_standardized <-
    endpoint_state$theta_tilde_standardized
  inference$correction_standardized <-
    endpoint_state$correction_standardized
  inference$pilot_residual_norm <- endpoint_state$pilot_residual_norm
  inference$score_second_moment <- endpoint_state$score_second_moment
  inference$psi_centered <- endpoint_state$psi_centered
  inference
}

.partial_max_inference <- function(
    endpoint_state,
    alpha,
    bootstrap_B,
    bootstrap_seed,
    bootstrap_multipliers,
    bootstrap_block_size) {
  inference <- max_test_hd(
    theta_tilde = endpoint_state$theta_tilde,
    target = numeric(length(endpoint_state$theta_tilde)),
    psi_centered = endpoint_state$psi_centered,
    alpha = alpha,
    bootstrap_B = bootstrap_B,
    bootstrap_seed = bootstrap_seed,
    multipliers = bootstrap_multipliers,
    block_size = bootstrap_block_size
  )
  inference$subtype <- "partial_debiased_diagonal_max"
  .attach_partial_endpoint_state(inference, endpoint_state)
}

.validate_threshold_eta <- function(threshold_eta, q) {
  threshold_eta <- as.numeric(threshold_eta)
  if (length(threshold_eta) != 1L || !is.finite(threshold_eta) ||
      threshold_eta <= 0 || threshold_eta >= q) {
    stop(
      "threshold_eta must be one finite value strictly between zero and q.",
      call. = FALSE
    )
  }
  threshold <- stats::qnorm(1 - threshold_eta / (2 * q))
  if (!is.finite(threshold) || threshold < 0) {
    stop("threshold_eta produced an invalid Gaussian threshold.", call. = FALSE)
  }
  list(eta = threshold_eta, threshold = threshold)
}

.default_threshold_eta <- function(q) {
  q <- as.integer(q)
  if (length(q) != 1L || is.na(q) || q < 1L) {
    stop("q must be one positive integer.", call. = FALSE)
  }
  min(5, q / 2)
}

.vanishing_threshold_eta <- function(q) {
  q <- as.integer(q)
  if (length(q) != 1L || is.na(q) || q < 1L) {
    stop("q must be one positive integer.", call. = FALSE)
  }
  if (q == 1L) 0.5 else 1 / log(q)
}

.threshold_partial_endpoint <- function(
    theta_tested,
    standard_error,
    prep,
    geometry,
    core_set,
    tested_set,
    threshold_eta) {
  n <- nrow(prep$X)
  p <- ncol(prep$X)
  q <- length(tested_set)
  theta_tested <- as.numeric(theta_tested)
  standard_error <- as.numeric(standard_error)
  if (length(theta_tested) != q || length(standard_error) != q ||
      any(!is.finite(theta_tested)) || any(!is.finite(standard_error)) ||
      any(standard_error <= 0)) {
    stop("Invalid debiased endpoint inputs.", call. = FALSE)
  }
  threshold_state <- .validate_threshold_eta(threshold_eta, q)
  residual_df <- n - as.integer(geometry$rank) - 1L
  if (length(residual_df) != 1L || is.na(residual_df) || residual_df <= 1L) {
    stop(
      "The residual degrees of freedom must exceed one for thresholding.",
      call. = FALSE
    )
  }
  studentization_correction <- sqrt(residual_df / n)
  studentized <- sqrt(residual_df) * theta_tested / standard_error
  retained <- abs(studentized) > threshold_state$threshold

  beta_standardized <- numeric(p)
  beta_standardized[tested_set] <-
    theta_tested * prep$x_scale[tested_set] * retained
  beta_standardized[core_set] <- as.numeric(
    geometry$coefficient_operator %*%
      (prep$y - prep$X[, tested_set, drop = FALSE] %*%
         beta_standardized[tested_set])
  )
  beta <- beta_standardized / prep$x_scale
  list(
    beta = beta,
    beta_standardized = beta_standardized,
    beta_tested_debiased = theta_tested,
    studentized = studentized,
    retained = retained,
    support_tested = tested_set[retained],
    eta = threshold_state$eta,
    threshold = threshold_state$threshold,
    expected_null_exceedances = threshold_state$eta,
    residual_df = residual_df,
    studentization_reference_n = n,
    studentization_correction = studentization_correction
  )
}

#' High-dimensional square-root-LASSO shrinkage toward an exact coordinate null
#'
#' This is a non-Ridge framework for a sparse linear model with `p > n`.
#' The columns are split in advance into a low-dimensional core `A` and a
#' high-dimensional tested block `B`. The recommended all-X endpoint centers
#' and scales both `X` and `y`, fixes the standardized intercept at zero, and
#' estimates every coefficient jointly in one square-root-LASSO problem.
#' Legacy partial endpoints remain available for reproducibility.
#' The submodel endpoint is the exact-null fit under `beta[B] = 0`; it is not
#' an oracle fit and is not obtained by truncating the full estimate.
#'
#' @param X Numeric `n` by `p` design matrix.
#' @param y Numeric response vector.
#' @param core_set Prespecified one-based indices of the low-dimensional core.
#'   It must contain fewer than `n` columns and must not be chosen from `y`
#'   without sample splitting.
#' @param tested_set Prespecified high-dimensional null block. By default it is
#'   the complement of `core_set`; together the two sets must partition `X`.
#' @param lambda Positive square-root-LASSO penalty. The default is
#'   `1.1 * sqrt(2 * log(2 * p) / n)` for the all-X endpoint and
#'   `1.1 * sqrt(2 * log(2 * q) / n)` for a legacy partial endpoint.
#' @param penalty_loadings Optional positive loadings for the tested block.
#'   The default is the empirical norm of each tested column after projection
#'   off the core.
#' @param all_x_penalty_loadings Optional length-`p` positive loadings for
#'   `full_endpoint = "all_x_square_root_lasso"`. The default penalizes all
#'   columns equally. No core coefficient is exempted from the FM penalty.
#' @param standardize_y Whether to center and scale `y` before fitting. It is
#'   fixed at `TRUE` for the all-X endpoint and must remain `FALSE` for the
#'   legacy endpoints. `X` is always centered and scaled. The standardized
#'   intercept is exactly zero; returned coefficients and prediction
#'   intercepts are back-transformed to the original data scale.
#' @param threshold_eta Positive marginal exceedance budget for the thresholded
#'   debiased full endpoint. The threshold is
#'   `qnorm(1 - threshold_eta / (2 * q))`. The default is `min(5, q / 2)`,
#'   which gives the prespecified finite-sample budget `5` in the advertised
#'   ultra-high-dimensional designs. This is a PFER-style tuning
#'   interpretation, not an FDR guarantee. The vanishing sequence
#'   `1 / log(q)` is retained for theoretical and sensitivity analyses. This
#'   argument must remain `NULL` for either unthresholded square-root-LASSO
#'   endpoint because no coefficient threshold is applied there.
#' @param assume_independent Must be `TRUE` for the default thresholded
#'   partial-debiased FM endpoint. It records the Gaussian-iid (or otherwise
#'   prespecified diagonal population precision) condition used by that
#'   endpoint. It is not required by either unthresholded square-root-LASSO
#'   endpoint or by the conditional-Gaussian exact-null test.
#' @param alpha Test level.
#' @param bootstrap_B Number of Gaussian multiplier draws.
#' @param bootstrap_seed Optional seed that leaves the caller RNG unchanged.
#' @param bootstrap_multipliers Optional reusable `n` by `bootstrap_B` matrix.
#'   Under the default conditional-Gaussian calibration it is accepted as a
#'   backward-compatible alias for `gaussian_draws`; under
#'   `test_calibration = "residual_multiplier"` it contains the multiplier
#'   draws used by the legacy test.
#' @param bootstrap_block_size Number of draws processed per C++ block.
#' @param max_iter Maximum primal-dual iterations.
#' @param tol Solver convergence tolerance.
#' @param check_every Iteration interval for diagnostics.
#' @param verbose Print periodic solver diagnostics.
#' @param full_beta_init Optional original-scale length-`p` warm start.
#' @param full_operator_norm Optional cached operator norm of the full
#'   standardized design for the all-X endpoint or of the residualized tested
#'   design for a legacy partial endpoint.
#' @param null_geometry Optional cached geometry returned in a previous fit on
#'   exactly the same design and partition.
#' @param rank_tolerance Relative numerical-rank tolerance for the core SVD.
#' @param epsilon Positive numerical floor in the Stein denominator.
#' @param test_calibration Test calibration. The default
#'   `"conditional_gaussian"` is the finite-sample exact-null Gaussian
#'   max-partial-t calibration. `"residual_multiplier"` retains the legacy
#'   partial-debiased multiplier maximum.
#' @param shrinkage_calibration Max-Stein calibration. The default
#'   `"inverse_moment"` uses the reciprocal bootstrap inverse second moment.
#'   `"second_moment"` retains the legacy bootstrap second-moment rule.
#' @param gaussian_draws Optional reusable `n` by `bootstrap_B` matrix for the
#'   conditional-Gaussian exact-null calibration.
#' @param full_endpoint Definition of the full-model endpoint. The default
#'   `"thresholded_partial_debiased"` preserves the Gaussian-iid framework.
#'   `"profiled_square_root_lasso"` uses the reprofiled partial
#'   square-root-LASSO pilot directly and is available for correlated fixed
#'   designs. `"all_x_square_root_lasso"` estimates all `p` coefficients in
#'   one centered-and-scaled square-root-LASSO problem and is the recommended
#'   endpoint when FM must use the complete design without an A/B-specific
#'   estimation step. Both unthresholded endpoints require
#'   `test_calibration = "conditional_gaussian"`.
#'
#' @return An object of class `hd_shrinkage_fit`. `FM` is the selected full
#'   endpoint, `SM` is the exact-null endpoint, and `PT`, `S`, `PS`, and `PPS`
#'   are test-guided estimators between those endpoints. `PPS` equals `FM`
#'   after rejection and otherwise equals `PS`; it is reported as a separate
#'   protected estimator and does not replace `PS`. `full_regularized` equals
#'   `FM` for either unthresholded square-root-LASSO endpoint.
#'
#' @details
#' For the recommended all-X endpoint, write `X*` and `y*` for columnwise and
#' responsewise centered-and-RMS-scaled training data. It solves once
#' \deqn{\widehat\beta^{FM,*}=\arg\min_{b\in\mathbb R^p}
#' \|y^*-X^*b\|_2/\sqrt n+\lambda\sum_{j=1}^p\omega_j|b_j|.}
#' Its standardized intercept is exactly zero. Coefficients and the prediction
#' intercept are then mapped back to the original scale. The A/B partition is
#' used only to construct the exact-null SM and the restriction test, so
#' changing that partition does not change the all-X FM for fixed `X`, `y`,
#' penalty, and solver settings.
#'
#' For a legacy partial endpoint, let `M_A` project orthogonally off the
#' core, `r_A = M_A y`, and `V = M_A X_B`. The pilot solves
#' \deqn{\hat b_B^{FM}=\arg\min_b
#' \|r_A-Vb\|_2/\sqrt n+\lambda\sum_j\omega_j|b_j|,}
#' and produces one-step coordinates `tilde beta_j` and standard errors. Let
#' `d = n - rank(X_A) - 1` be the residual degrees of freedom after centering
#' and projecting off the core. With
#' \deqn{Z_j=\sqrt d\,\widetilde\beta_j/\widehat{se}_j,\qquad
#' \tau_\eta=\Phi^{-1}(1-\eta/(2q)),}
#' the full endpoint is
#' \deqn{\widehat\beta_{B,j}^{FM}=\widetilde\beta_j
#' 1\{|Z_j|>\tau_\eta\},\qquad
#' \widehat\beta_A^{FM}=X_A^+(y-X_B\widehat\beta_B^{FM}).}
#' Under `full_endpoint = "profiled_square_root_lasso"`, the tested block is
#' instead the optimizer `hat b_B` in the first display and
#' \deqn{\widehat\beta_A^{FM}=X_A^+(y-X_B\hat b_B).}
#' This alternative does not debias or threshold the tested coefficients and
#' therefore does not invoke a diagonal precision approximation.
#' The exact-null endpoint is
#' \deqn{\hat b^{SM}=(X_A^+y,0_q).}
#' The Moore--Penrose operator is computed from a checked economy SVD; no OLS
#' Gram inverse and no `p` by `p` inverse is used.
#'
#' Conditional on the fixed design, the default test starts from the exact-null
#' residual, applies the same core projection to every Gaussian draw, and
#' transforms the maximum normalized score monotonically to the maximum
#' partial-t scale. Its Monte Carlo p-value is therefore unchanged by the
#' transformation. No direct chi-square limit and no `q` by `q` covariance
#' inverse is asserted. By default the max-Stein factor is
#' `1 - kappa_star / T_max^2`, where
#' `kappa_star = 1 / mean(1 / (T_max_star)^2)`; the positive-part factor
#' truncates this at zero. With `R` denoting the rejection indicator, the
#' protected weight is `R + (1 - R) * positive_weight`. Classical
#' James--Stein dominance is not claimed for either estimator. Finite-sample
#' conditional exactness requires homoskedastic Gaussian errors and a
#' prespecified core, but it does not require independent design columns.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(50 * 100), 50, 100)
#' beta <- c(1.2, -1, 0.8, rep(0, 97))
#' y <- drop(X %*% beta + rnorm(50))
#' fit <- fit_partial_sqrt_lasso_shrinkage(
#'   X, y, core_set = 1:3, tested_set = 4:100,
#'   assume_independent = TRUE, bootstrap_B = 99,
#'   max_iter = 5000, tol = 1e-4
#' )
#' coef(fit, method = "PS")
#' coef(fit, method = "PPS")
#' fit_all_x <- fit_partial_sqrt_lasso_shrinkage(
#'   X, y, core_set = 1:3, tested_set = 4:100,
#'   bootstrap_B = 99, max_iter = 5000, tol = 1e-4,
#'   full_endpoint = "all_x_square_root_lasso"
#' )
#'
#' @export
fit_partial_sqrt_lasso_shrinkage <- function(
    X,
    y,
    core_set,
    tested_set = NULL,
    lambda = NULL,
    penalty_loadings = NULL,
    all_x_penalty_loadings = NULL,
    standardize_y = NULL,
    threshold_eta = NULL,
    assume_independent = FALSE,
    alpha = 0.05,
    bootstrap_B = 999L,
    bootstrap_seed = NULL,
    bootstrap_multipliers = NULL,
    bootstrap_block_size = 256L,
    max_iter = 30000L,
    tol = 5e-6,
    check_every = 50L,
    verbose = FALSE,
    full_beta_init = NULL,
    full_operator_norm = NULL,
    null_geometry = NULL,
    rank_tolerance = 1e-10,
    epsilon = 1e-10,
    test_calibration = c("conditional_gaussian", "residual_multiplier"),
    shrinkage_calibration = c("inverse_moment", "second_moment"),
    gaussian_draws = NULL,
    full_endpoint = c(
      "thresholded_partial_debiased", "profiled_square_root_lasso",
      "all_x_square_root_lasso"
    )) {
  call <- match.call()
  full_endpoint <- match.arg(full_endpoint)
  all_x_endpoint <- identical(full_endpoint, "all_x_square_root_lasso")
  standardize_y <- standardize_y %||% all_x_endpoint
  if (!is.logical(standardize_y) || length(standardize_y) != 1L ||
      is.na(standardize_y)) {
    stop("standardize_y must be TRUE or FALSE.", call. = FALSE)
  }
  if (all_x_endpoint && !isTRUE(standardize_y)) {
    stop(
      "standardize_y must be TRUE for the all-X full endpoint so X and y ",
      "share the frozen centered/RMS-scaled, zero-intercept convention.",
      call. = FALSE
    )
  }
  if (isTRUE(standardize_y) && !all_x_endpoint) {
    stop(
      "standardize_y = TRUE is currently reserved for the all-X full ",
      "endpoint so legacy endpoint scales remain unchanged.",
      call. = FALSE
    )
  }
  prep <- .prepare_design(X, y, scale_y = isTRUE(standardize_y))
  n <- nrow(prep$X)
  p <- ncol(prep$X)
  partition <- .partial_partition(core_set, tested_set, p)
  core_set <- partition$core_set
  tested_set <- partition$tested_set
  p1 <- length(core_set)
  q <- length(tested_set)
  if (p1 >= n) {
    stop("length(core_set) must be strictly smaller than n.", call. = FALSE)
  }
  test_calibration <- match.arg(test_calibration)
  shrinkage_calibration <- match.arg(shrinkage_calibration)
  thresholded_debiased_endpoint <- identical(
    full_endpoint, "thresholded_partial_debiased"
  )
  profiled_endpoint <- identical(
    full_endpoint, "profiled_square_root_lasso"
  )
  if (thresholded_debiased_endpoint && !isTRUE(assume_independent)) {
    stop(
      "Set assume_independent = TRUE only when diagonal population precision ",
      "is prespecified, as in the Gaussian-iid design.",
      call. = FALSE
    )
  }
  if (!thresholded_debiased_endpoint &&
      !identical(test_calibration, "conditional_gaussian")) {
    stop(
      "The selected non-debiased full endpoint requires ",
      "test_calibration = \"conditional_gaussian\".",
      call. = FALSE
    )
  }
  if (!thresholded_debiased_endpoint && !is.null(threshold_eta)) {
    stop(
      "threshold_eta is not used by the selected square-root-LASSO ",
      "full endpoint; leave it NULL.",
      call. = FALSE
    )
  }
  lambda_dimension <- if (all_x_endpoint) p else q
  lambda <- lambda %||%
    (1.1 * sqrt(2 * log(2 * lambda_dimension) / n))
  lambda <- .check_positive_scalar(lambda, "lambda")
  if (thresholded_debiased_endpoint) {
    threshold_eta <- threshold_eta %||% .default_threshold_eta(q)
    threshold_state <- .validate_threshold_eta(threshold_eta, q)
    threshold_eta <- threshold_state$eta
  }
  rank_tolerance <- .check_positive_scalar(
    rank_tolerance, "rank_tolerance"
  )
  epsilon <- .check_positive_scalar(epsilon, "epsilon")
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie strictly between zero and one.", call. = FALSE)
  }

  if (is.null(null_geometry)) {
    null_geometry <- .make_partial_null_geometry(
      X = prep$X,
      core_set = core_set,
      tested_set = tested_set,
      rank_tolerance = rank_tolerance
    )
  } else {
    .validate_partial_null_geometry(
      null_geometry, prep$X, core_set, tested_set
    )
  }
  null_state <- cpp_partial_null_apply(
    y = prep$y,
    coefficient_operator = null_geometry$coefficient_operator,
    orthonormal_basis = null_geometry$orthonormal_basis,
    residualized_tested = null_geometry$residualized_tested
  )

  if (all_x_endpoint && !is.null(penalty_loadings)) {
    stop(
      "penalty_loadings is specific to the legacy partial endpoints; use ",
      "all_x_penalty_loadings for the all-X endpoint.",
      call. = FALSE
    )
  }
  if (!all_x_endpoint && !is.null(all_x_penalty_loadings)) {
    stop(
      "all_x_penalty_loadings is used only by full_endpoint = ",
      "\"all_x_square_root_lasso\".",
      call. = FALSE
    )
  }
  if (all_x_endpoint) {
    all_x_penalty_loadings <- all_x_penalty_loadings %||% rep(1, p)
    all_x_penalty_loadings <- as.numeric(all_x_penalty_loadings)
    if (length(all_x_penalty_loadings) != p ||
        any(!is.finite(all_x_penalty_loadings)) ||
        any(all_x_penalty_loadings <= 0)) {
      stop(
        "all_x_penalty_loadings must contain p finite positive values.",
        call. = FALSE
      )
    }
    penalty_loadings <- NULL
  } else {
    penalty_loadings <- penalty_loadings %||%
      sqrt(null_geometry$residualized_second_moment)
    penalty_loadings <- as.numeric(penalty_loadings)
    if (length(penalty_loadings) != q || any(!is.finite(penalty_loadings)) ||
        any(penalty_loadings <= 0)) {
      stop(
        "penalty_loadings must contain q finite positive values.",
        call. = FALSE
      )
    }
  }
  beta_init_standardized <- numeric(p)
  if (!is.null(full_beta_init)) {
    full_beta_init <- as.numeric(full_beta_init)
    if (length(full_beta_init) != p || any(!is.finite(full_beta_init))) {
      stop("full_beta_init must contain p finite values.", call. = FALSE)
    }
    beta_init_standardized <-
      full_beta_init * prep$x_scale / prep$y_scale
  }

  if (all_x_endpoint) {
    endpoint_solver <- .sqrt_lasso_standardized(
      X = prep$X,
      y = prep$y,
      lambda = lambda,
      penalty_factor = all_x_penalty_loadings,
      beta_init = beta_init_standardized,
      max_iter = max_iter,
      tol = tol,
      check_every = check_every,
      verbose = verbose,
      operator_norm_value = full_operator_norm
    )
    profiled_solver <- NULL
    beta_pilot_standardized <- as.numeric(endpoint_solver$beta)
  } else {
    beta_init_tested <- beta_init_standardized[tested_set]
    profiled_solver <- .sqrt_lasso_standardized(
      X = null_geometry$residualized_tested,
      y = null_state$residual,
      lambda = lambda,
      penalty_factor = penalty_loadings,
      beta_init = beta_init_tested,
      max_iter = max_iter,
      tol = tol,
      check_every = check_every,
      verbose = verbose,
      operator_norm_value = full_operator_norm
    )
    endpoint_solver <- profiled_solver
    beta_pilot_standardized <- numeric(p)
    beta_pilot_standardized[tested_set] <- profiled_solver$beta
    beta_pilot_standardized[core_set] <- as.numeric(
      null_geometry$coefficient_operator %*%
        (prep$y - prep$X[, tested_set, drop = FALSE] %*%
           profiled_solver$beta)
    )
  }
  beta_submodel_standardized <- numeric(p)
  beta_submodel_standardized[core_set] <- null_state$beta_core

  if (!isTRUE(endpoint_solver$converged)) {
    warning(
      "The selected square-root LASSO reached max_iter before ",
      "meeting the requested tolerance."
    )
  }

  endpoint_state <- NULL
  if (thresholded_debiased_endpoint) {
    endpoint_state <- .partial_debiased_endpoint_state(
      null_state = null_state,
      geometry = null_geometry,
      beta_tested_standardized = profiled_solver$beta,
      tested_scale = prep$x_scale[tested_set]
    )
  }
  if (identical(test_calibration, "conditional_gaussian")) {
    if (!is.null(gaussian_draws) && !is.null(bootstrap_multipliers)) {
      stop(
        "Supply only one of gaussian_draws and bootstrap_multipliers under ",
        "conditional-Gaussian calibration.",
        call. = FALSE
      )
    }
    conditional_draws <- gaussian_draws %||% bootstrap_multipliers
    inference <- .gaussian_score_max_inference(
      null_state = null_state,
      geometry = null_geometry,
      alpha = alpha,
      bootstrap_B = bootstrap_B,
      bootstrap_seed = bootstrap_seed,
      gaussian_draws = conditional_draws,
      bootstrap_block_size = bootstrap_block_size
    )
    inference <- .as_max_partial_t_inference(inference)
    if (thresholded_debiased_endpoint) {
      inference <- .attach_partial_endpoint_state(inference, endpoint_state)
    }
    inference$draw_source <- if (is.null(gaussian_draws)) {
      if (is.null(bootstrap_multipliers)) "generated_gaussian" else
        "bootstrap_multipliers_alias"
    } else {
      "gaussian_draws"
    }
  } else {
    if (!is.null(gaussian_draws)) {
      stop(
        "gaussian_draws is only used by conditional-Gaussian calibration.",
        call. = FALSE
      )
    }
    inference <- .partial_max_inference(
      endpoint_state = endpoint_state,
      alpha = alpha,
      bootstrap_B = bootstrap_B,
      bootstrap_seed = bootstrap_seed,
      bootstrap_multipliers = bootstrap_multipliers,
      bootstrap_block_size = bootstrap_block_size
    )
  }
  inference <- .set_max_shrinkage_calibration(
    inference, calibration_type = shrinkage_calibration
  )
  inference$test_calibration <- test_calibration
  coefficient_backscale <- prep$y_scale / prep$x_scale
  beta_pilot <- beta_pilot_standardized * coefficient_backscale
  thresholded_endpoint <- NULL
  beta_debiased_standardized <- NULL
  beta_debiased <- NULL
  if (thresholded_debiased_endpoint) {
    thresholded_endpoint <- .threshold_partial_endpoint(
      theta_tested = inference$theta_tilde,
      standard_error = inference$standard_error,
      prep = prep,
      geometry = null_geometry,
      core_set = core_set,
      tested_set = tested_set,
      threshold_eta = threshold_eta
    )
    beta_full <- thresholded_endpoint$beta
    beta_full_standardized <- thresholded_endpoint$beta_standardized

    beta_debiased_standardized <- numeric(p)
    beta_debiased_standardized[tested_set] <-
      inference$theta_tilde * prep$x_scale[tested_set]
    beta_debiased_standardized[core_set] <- as.numeric(
      null_geometry$coefficient_operator %*%
        (prep$y - prep$X[, tested_set, drop = FALSE] %*%
           beta_debiased_standardized[tested_set])
    )
    beta_debiased <- beta_debiased_standardized / prep$x_scale
  } else {
    beta_full <- beta_pilot
    beta_full_standardized <- beta_pilot_standardized
  }
  beta_submodel <- beta_submodel_standardized * coefficient_backscale
  shrinkage <- .make_max_shrinkage(
    beta_full = beta_full,
    beta_restricted = beta_submodel,
    statistic = inference$statistic,
    reject = inference$reject,
    epsilon = epsilon,
    shrinkage_calibration = inference$shrinkage_calibration,
    shrinkage_calibration_type =
      inference$shrinkage_calibration_type
  )

  all_beta <- list(
    full = beta_full,
    submodel = beta_submodel,
    restricted = beta_submodel,
    full_regularized = beta_pilot,
    submodel_regularized = beta_submodel,
    full_debiased = beta_debiased,
    submodel_debiased = beta_submodel,
    preliminary_test = shrinkage$preliminary_test,
    stein = shrinkage$stein,
    positive_part = shrinkage$positive_part,
    protected_positive_part = shrinkage$protected_positive_part
  )
  intercept <- lapply(all_beta, function(beta) {
    if (is.null(beta)) return(NULL)
    prep$y_center - sum(prep$x_center * beta)
  })
  standardized_intercept <- lapply(all_beta, function(beta) {
    if (is.null(beta)) return(NULL)
    0
  })

  full_solver <- endpoint_solver
  full_solver$beta_profiled <- if (!all_x_endpoint) {
    profiled_solver$beta
  } else {
    NULL
  }
  full_solver$beta <- beta_full_standardized
  full_solver$beta_pilot <- beta_pilot_standardized
  full_solver$beta_core <- beta_full_standardized[core_set]
  full_solver$beta_tested <- beta_full_standardized[tested_set]
  full_solver$penalty_loadings <- if (all_x_endpoint) {
    all_x_penalty_loadings
  } else {
    penalty_loadings
  }
  full_solver$penalty_scope <- if (all_x_endpoint) "all_X" else "tested_block"
  full_solver$standardized_response <- isTRUE(standardize_y)
  full_solver$support_pilot_tested <- tested_set[
    abs(beta_pilot[tested_set]) > 1e-8
  ]
  if (thresholded_debiased_endpoint) {
    full_solver$beta_debiased <- beta_debiased_standardized
    full_solver$beta_tested_debiased <-
      thresholded_endpoint$beta_tested_debiased
    full_solver$studentized <- thresholded_endpoint$studentized
    full_solver$threshold_eta <- thresholded_endpoint$eta
    full_solver$threshold <- thresholded_endpoint$threshold
    full_solver$expected_null_exceedances <-
      thresholded_endpoint$expected_null_exceedances
    full_solver$threshold_residual_df <- thresholded_endpoint$residual_df
    full_solver$threshold_studentization_reference_n <-
      thresholded_endpoint$studentization_reference_n
    full_solver$threshold_studentization_correction <-
      thresholded_endpoint$studentization_correction
    full_solver$support_tested <- thresholded_endpoint$support_tested
    full_solver$construction <-
      "thresholded_partial_debiased_square_root_lasso"
  } else if (profiled_endpoint) {
    full_solver$support_tested <- full_solver$support_pilot_tested
    full_solver$construction <- "profiled_partial_square_root_lasso"
  } else {
    full_solver$support_tested <- full_solver$support_pilot_tested
    full_solver$support_all <- which(abs(beta_pilot) > 1e-8)
    full_solver$construction <- "all_x_square_root_lasso"
  }
  restricted_solver <- list(
    beta = beta_submodel_standardized,
    beta_core = beta_submodel_standardized[core_set],
    beta_tested = beta_submodel_standardized[tested_set],
    converged = TRUE,
    iterations = 1L,
    objective = null_state$residual_norm / sqrt(n),
    restriction_violation = max(abs(beta_submodel[tested_set])),
    construction = "exact_null_refit_svd",
    rank = null_geometry$rank,
    condition_number = null_geometry$condition_number,
    residual_norm = null_state$residual_norm
  )
  restriction <- list(
    type = "coordinate",
    M = tested_set,
    core_set = core_set,
    tested_set = tested_set,
    C = NULL,
    t = numeric(q),
    zero_block = TRUE,
    construction = "exact_null_refit_svd",
    geometry = "profiled_core_SVD",
    exact_endpoint_restriction = TRUE
  )
  geometry_diagnostics <- list(
    rank = null_geometry$rank,
    singular_values = null_geometry$singular_values,
    condition_number = null_geometry$condition_number,
    orthogonality_error = null_geometry$orthogonality_error,
    rank_tolerance = null_geometry$rank_tolerance
  )
  if (thresholded_debiased_endpoint) {
    precision_output <- list(
      method = "prespecified_diagonal_after_core_residualization",
      assume_independent = TRUE
    )
    threshold_output <- list(
      type = "marginal_gaussian_pfer",
      eta = thresholded_endpoint$eta,
      value = thresholded_endpoint$threshold,
      expected_null_exceedances =
        thresholded_endpoint$expected_null_exceedances,
      residual_df = thresholded_endpoint$residual_df,
      studentization_reference_n =
        thresholded_endpoint$studentization_reference_n,
      studentization_correction =
        thresholded_endpoint$studentization_correction,
      studentization = "residual_df_adjusted_partial_debiased",
      studentized = thresholded_endpoint$studentized,
      retained = thresholded_endpoint$retained,
      support_tested = thresholded_endpoint$support_tested
    )
  } else {
    precision_output <- list(
      method = if (all_x_endpoint) {
        "not_required_for_all_x_square_root_lasso_endpoint"
      } else {
        "not_required_for_profiled_square_root_lasso_endpoint"
      },
      assume_independent = isTRUE(assume_independent),
      required = FALSE
    )
    threshold_output <- list(
      type = "not_applicable",
      eta = NA_real_,
      value = NA_real_,
      expected_null_exceedances = NA_real_,
      residual_df = inference$residual_df,
      studentization_reference_n = n,
      studentization_correction = NA_real_,
      studentization = "not_applied",
      studentized = NULL,
      retained = NULL,
      support_tested = full_solver$support_tested
    )
  }

  output <- list(
    call = call,
    dimensions = c(n = n, p = p, p1 = p1, p2 = q, m = q, q = q),
    beta = all_beta,
    intercept = intercept,
    standardized_intercept = standardized_intercept,
    endpoint = "partial_sqrt_lasso_null",
    test = if (identical(test_calibration, "conditional_gaussian")) {
      "max_partial_t"
    } else {
      "max"
    },
    test_calibration = test_calibration,
    reject = inference$reject,
    alpha = alpha,
    inference = inference,
    shrinkage = shrinkage,
    full_solver = full_solver,
    restricted_solver = restricted_solver,
    restriction_projection = restricted_solver,
    restriction = restriction,
    precision = precision_output,
    endpoint_precision = NULL,
    threshold = threshold_output,
    score_representation = if (
      identical(test_calibration, "conditional_gaussian")
    ) {
      "conditional_gaussian_exact_null_partial_t"
    } else {
      "partial_debiased_diagonal"
    },
    preprocessing = prep[c(
      "x_center", "x_scale", "y_center", "y_scale",
      "standardized_intercept"
    )],
    standardized_restriction = list(
      tested_set = tested_set,
      target = numeric(q)
    ),
    null_geometry = null_geometry,
    geometry_diagnostics = geometry_diagnostics
  )
  if (!thresholded_debiased_endpoint) {
    output$full_endpoint <- full_endpoint
  }
  class(output) <- "hd_shrinkage_fit"
  output
}

#' Conditional-Gaussian exact-null maximum partial-t test
#'
#' This convenience function runs exactly the same full-model pilot and
#' maximum test used by
#' [fit_partial_sqrt_lasso_shrinkage()] and returns the inference component.
#' The default is the conditional-Gaussian exact-null maximum partial-t test;
#' it is not a Wald, chi-square, or joint partial-F test.
#'
#' @inheritParams fit_partial_sqrt_lasso_shrinkage
#'
#' @return A maximum-test list with the fitted pilot attached as
#'   `full_model_beta` and `full_solver`.
#' @export
max_score_test_hd <- function(
    X,
    y,
    core_set,
    tested_set = NULL,
    lambda = NULL,
    penalty_loadings = NULL,
    all_x_penalty_loadings = NULL,
    standardize_y = NULL,
    threshold_eta = NULL,
    assume_independent = FALSE,
    alpha = 0.05,
    bootstrap_B = 999L,
    bootstrap_seed = NULL,
    bootstrap_multipliers = NULL,
    bootstrap_block_size = 256L,
    max_iter = 30000L,
    tol = 5e-6,
    check_every = 50L,
    verbose = FALSE,
    full_beta_init = NULL,
    full_operator_norm = NULL,
    null_geometry = NULL,
    rank_tolerance = 1e-10,
    test_calibration = c("conditional_gaussian", "residual_multiplier"),
    shrinkage_calibration = c("inverse_moment", "second_moment"),
    gaussian_draws = NULL,
    full_endpoint = c(
      "thresholded_partial_debiased", "profiled_square_root_lasso",
      "all_x_square_root_lasso"
    )) {
  fit <- fit_partial_sqrt_lasso_shrinkage(
    X = X,
    y = y,
    core_set = core_set,
    tested_set = tested_set,
    lambda = lambda,
    penalty_loadings = penalty_loadings,
    all_x_penalty_loadings = all_x_penalty_loadings,
    standardize_y = standardize_y,
    threshold_eta = threshold_eta,
    assume_independent = assume_independent,
    alpha = alpha,
    bootstrap_B = bootstrap_B,
    bootstrap_seed = bootstrap_seed,
    bootstrap_multipliers = bootstrap_multipliers,
    bootstrap_block_size = bootstrap_block_size,
    max_iter = max_iter,
    tol = tol,
    check_every = check_every,
    verbose = verbose,
    full_beta_init = full_beta_init,
    full_operator_norm = full_operator_norm,
    null_geometry = null_geometry,
    rank_tolerance = rank_tolerance,
    test_calibration = test_calibration,
    shrinkage_calibration = shrinkage_calibration,
    gaussian_draws = gaussian_draws,
    full_endpoint = full_endpoint
  )
  output <- fit$inference
  output$full_model_beta <- fit$beta$full
  output$full_solver <- fit$full_solver
  output$core_set <- fit$restriction$core_set
  output$tested_set <- fit$restriction$tested_set
  if (!is.null(fit$full_endpoint)) {
    output$full_endpoint <- fit$full_endpoint
  }
  class(output) <- c("hd_partial_max_test", "list")
  output
}
