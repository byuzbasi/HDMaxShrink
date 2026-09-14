#' Restriction-adaptive shrinkage in a high-dimensional linear model
#'
#' `fit_hd_shrinkage()` fits one full square-root lasso with all columns of
#' `X`. The submodel endpoint is the restriction projection of that full
#' estimator; it is never re-fitted on a reduced design. The endpoints can be
#' the regularized full/projection pair or, for a zero-coordinate restriction,
#' a full one-step debiased estimator and its projection. A debiased Wald
#' statistic is available for fixed, small restriction dimension. A
#' multiplier-bootstrap maximum test is used for a large coordinate block and
#' does not form a `q` by `q` covariance matrix. Neither branch uses ordinary,
#' ridge, or generalized ridge regression.
#'
#' @param X Numeric `n` by `p` design matrix.
#' @param y Numeric response vector of length `n`.
#' @param restricted_set One-based indices `M` involved in the restriction.
#' @param C Optional `q` by `length(M)` contrast matrix. `NULL` means the
#'   coordinate restriction `beta[M] = t`, whose projection is implemented
#'   without constructing a `q` by `q` matrix.
#' @param t Null target. It defaults to zero.
#' @param endpoint Shrinkage endpoints. `"regularized"` uses the full
#'   square-root lasso and its restriction projection. `"debiased_models"`
#'   computes one full-model debiased estimator and projects it to the
#'   coordinate zero subspace. The latter currently requires a coordinate
#'   zero-block restriction.
#' @param lambda Square-root lasso penalty.
#' @param lambda_node Nodewise lasso penalty.
#' @param test One of `"auto"`, `"wald"`, or `"max"`. Automatic selection
#'   uses Wald only when `q <= fixed_q_max`.
#' @param precision One of `"nodewise"`, `"diagonal"`, or `"user"`.
#' @param theta_rows User precision rows for `precision = "user"`.
#' @param precision_fit Optional reusable result from
#'   [estimate_precision_rows()].
#' @param full_model_precision_fit Optional reusable all-row precision result
#'   for the standardized full design when `endpoint = "debiased_models"`.
#' @param submodel_precision_fit Deprecated compatibility argument. It is
#'   ignored because no separate submodel debiasing fit is performed.
#' @param assume_independent Required for diagonal precision. This is suitable
#'   for a design whose diagonal population precision is assumed in advance,
#'   such as the Gaussian iid simulation design.
#' @param alpha Test level.
#' @param bootstrap_B Number of multiplier draws for the maximum test.
#' @param bootstrap_seed Optional seed used without changing the caller's RNG
#'   state.
#' @param bootstrap_multipliers Optional reusable `n` by `bootstrap_B`
#'   Gaussian multiplier matrix.
#' @param bootstrap_block_size Number of multiplier draws processed in one C++
#'   block.
#' @param fixed_q_max Largest restriction dimension selected for Wald by
#'   `test = "auto"`.
#' @param penalty_factor Nonnegative square-root lasso penalty factors.
#' @param max_iter Maximum primal-dual iterations.
#' @param tol Solver tolerance.
#' @param check_every Iteration interval for convergence diagnostics.
#' @param verbose Print periodic solver diagnostics.
#' @param full_beta_init Optional original-scale warm start for the full fit.
#' @param restricted_beta_init Deprecated compatibility argument. It is
#'   ignored because the restricted endpoint is projected from the full fit.
#' @param full_operator_norm Optional reusable full-design operator norm.
#' @param restricted_operator_norm Deprecated compatibility argument. It is
#'   ignored because no restricted-design solver is run.
#' @param epsilon Positive numerical floor used only in a shrinkage denominator.
#'
#' @return An object of class `hd_shrinkage_fit` containing the regularized and
#'   selected full/submodel endpoints, the available adaptive estimators, full
#'   solver and projection diagnostics, precision directions, and inference
#'   output.
#'
#' @details
#' Let `A beta = t` denote the restriction and let `V` denote the projection
#' geometry. The conceptual restricted endpoint is
#' \deqn{\hat\beta^{SM}=\hat\beta^F-
#' V A^T(A V A^T)^\dagger(A\hat\beta^F-t).}
#' The implementation uses Euclidean geometry after internal column
#' standardization. For `beta[M] = t`, this simply replaces the selected
#' coordinates of the all-`X` full estimate by `t`; all other coordinates are
#' inherited exactly from the full fit.
#'
#' For a large coordinate restriction, the max-calibrated estimators are
#' \deqn{\hat\beta^{S_\infty}=\hat\beta^{SM}+
#' (1-a^*/(T_{\max}^2\vee\epsilon))
#' (\hat\beta^{FM}-\hat\beta^{SM})}
#' and its positive-part version. The separately reported protected version
#' equals the full endpoint when the max test rejects and otherwise equals the
#' positive-part estimator. Here `a*` is the conditional mean of the
#' squared multiplier-bootstrap maximum. These estimators are experimental;
#' the package does not assert classical James--Stein dominance for them.
#'
#' The restricted set must be specified independently of the response used for
#' inference. If it is selected from data, sample splitting or external
#' information is required for the advertised test interpretation.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(40 * 70), 40, 70)
#' beta <- c(rep(1, 4), rep(0, 66))
#' y <- drop(X %*% beta + rnorm(40))
#' fit <- fit_hd_shrinkage(
#'   X, y,
#'   restricted_set = 5:70,
#'   test = "max",
#'   precision = "diagonal",
#'   assume_independent = TRUE,
#'   bootstrap_B = 99,
#'   max_iter = 5000,
#'   tol = 1e-4
#' )
#' coef(fit, method = "PS")
#' coef(fit, method = "PPS")
#'
#' @export
fit_hd_shrinkage <- function(
    X,
    y,
    restricted_set,
    C = NULL,
    t = NULL,
    endpoint = c("regularized", "debiased_models"),
    lambda = NULL,
    lambda_node = NULL,
    test = c("auto", "wald", "max"),
    precision = c("nodewise", "diagonal", "user"),
    theta_rows = NULL,
    precision_fit = NULL,
    full_model_precision_fit = NULL,
    submodel_precision_fit = NULL,
    assume_independent = FALSE,
    alpha = 0.05,
    bootstrap_B = 999L,
    bootstrap_seed = NULL,
    bootstrap_multipliers = NULL,
    bootstrap_block_size = 256L,
    fixed_q_max = 10L,
    penalty_factor = NULL,
    max_iter = 30000L,
    tol = 5e-6,
    check_every = 50L,
    verbose = FALSE,
    full_beta_init = NULL,
    restricted_beta_init = NULL,
    full_operator_norm = NULL,
    restricted_operator_norm = NULL,
    epsilon = 1e-10) {
  call <- match.call()
  test <- match.arg(test)
  precision <- match.arg(precision)
  endpoint <- match.arg(endpoint)
  prep <- .prepare_design(X, y)
  n <- nrow(prep$X)
  p <- ncol(prep$X)
  restricted_set <- .validate_indices(restricted_set, p)
  m <- length(restricted_set)

  coordinate_restriction <- is.null(C)
  if (coordinate_restriction) {
    q <- m
  } else {
    C <- .validate_numeric_matrix(C, "C")
    if (ncol(C) != m) {
      stop("ncol(C) must equal length(restricted_set).", call. = FALSE)
    }
    q <- nrow(C)
    if (qr(C)$rank != q) {
      stop("C must have full row rank.", call. = FALSE)
    }
  }
  t <- t %||% numeric(q)
  t <- as.numeric(t)
  if (length(t) != q || any(!is.finite(t))) {
    stop("t must contain q finite entries.", call. = FALSE)
  }
  zero_block <- coordinate_restriction && all(t == 0)
  if (endpoint == "debiased_models" && !zero_block) {
    stop(
      "endpoint = 'debiased_models' currently requires the coordinate ",
      "restriction beta[restricted_set] = 0.",
      call. = FALSE
    )
  }

  lambda <- lambda %||% (1.1 * sqrt(2 * log(2 * p) / n))
  .check_positive_scalar(lambda, "lambda")
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie strictly between zero and one.", call. = FALSE)
  }
  fixed_q_max <- as.integer(fixed_q_max)
  if (length(fixed_q_max) != 1L || is.na(fixed_q_max) || fixed_q_max < 1L) {
    stop("fixed_q_max must be a positive integer.", call. = FALSE)
  }
  epsilon <- .check_positive_scalar(epsilon, "epsilon")

  penalty_factor <- penalty_factor %||% rep(1, p)
  penalty_factor <- as.numeric(penalty_factor)
  if (length(penalty_factor) != p || any(!is.finite(penalty_factor)) ||
      any(penalty_factor < 0)) {
    stop(
      "penalty_factor must contain p finite nonnegative values.",
      call. = FALSE
    )
  }
  standardize_warm_start <- function(value, name) {
    if (is.null(value)) return(numeric(p))
    value <- as.numeric(value)
    if (length(value) != p || any(!is.finite(value))) {
      stop(name, " must contain p finite coefficients.", call. = FALSE)
    }
    value * prep$x_scale
  }
  full_init_standardized <- standardize_warm_start(
    full_beta_init, "full_beta_init"
  )

  full_solver <- .sqrt_lasso_standardized(
    X = prep$X,
    y = prep$y,
    lambda = lambda,
    penalty_factor = penalty_factor,
    beta_init = full_init_standardized,
    max_iter = max_iter,
    tol = tol,
    check_every = check_every,
    verbose = verbose,
    operator_norm_value = full_operator_norm
  )

  C_standardized <- NULL
  if (coordinate_restriction) {
    projection_target <- t * prep$x_scale[restricted_set]
  } else {
    C_standardized <- sweep(
      C, 2L, prep$x_scale[restricted_set], "/"
    )
    projection_target <- t
  }
  projection <- .project_full_standardized(
    beta_full = full_solver$beta,
    restricted_set = restricted_set,
    coordinate_restriction = coordinate_restriction,
    C_standardized = C_standardized,
    target = projection_target
  )
  standardized_restriction <- projection$standardized_restriction
  restricted_solver <- .make_projection_state(
    full_solver = full_solver,
    beta_projected = projection$beta,
    X = prep$X,
    y = prep$y,
    lambda = lambda,
    penalty_factor = penalty_factor,
    restriction_violation = projection$restriction_violation
  )

  if (!isTRUE(full_solver$converged)) {
    warning("Full square-root lasso reached max_iter before convergence.")
  }

  full_debias_state <- NULL
  beta_full_debiased <- NULL
  beta_submodel_debiased <- NULL
  endpoint_precision <- list(full = NULL, submodel = NULL)
  reuse_full_score <- FALSE

  if (endpoint == "debiased_models") {
    if (is.null(full_model_precision_fit)) {
      if (precision == "user") {
        stop(
          "precision = 'user' with debiased model endpoints requires ",
          "full_model_precision_fit.",
          call. = FALSE
        )
      }
      full_model_precision_fit <- estimate_precision_rows(
        X = prep$X,
        targets = seq_len(p),
        method = precision,
        lambda_node = lambda_node,
        assume_independent = assume_independent
      )
    }
    .validate_precision_fit_for_targets(
      full_model_precision_fit,
      targets = seq_len(p),
      p = p,
      name = "full_model_precision_fit"
    )
    if (!identical(full_model_precision_fit$method, precision)) {
      stop(
        "full_model_precision_fit must use the selected precision method.",
        call. = FALSE
      )
    }

    full_debias_state <- .debiased_coordinate_model(
      X = prep$X,
      y = prep$y,
      beta_standardized = full_solver$beta,
      x_scale = prep$x_scale,
      precision_fit = full_model_precision_fit
    )
    beta_full_debiased <- full_debias_state$theta_tilde
    debiased_projection <- .project_full_standardized(
      beta_full = beta_full_debiased * prep$x_scale,
      restricted_set = restricted_set,
      coordinate_restriction = coordinate_restriction,
      C_standardized = C_standardized,
      target = projection_target
    )
    beta_submodel_debiased <-
      debiased_projection$beta / prep$x_scale
    endpoint_precision <- list(
      full = full_model_precision_fit,
      submodel = full_model_precision_fit,
      shared_full_fit = TRUE
    )

    if (is.null(precision_fit)) {
      precision_fit <- .subset_precision_fit(
        full_model_precision_fit, restricted_set
      )
      reuse_full_score <- TRUE
    }
  }

  if (is.null(precision_fit)) {
    precision_fit <- estimate_precision_rows(
      X = prep$X,
      targets = restricted_set,
      method = precision,
      lambda_node = lambda_node,
      theta_rows = theta_rows,
      assume_independent = assume_independent
    )
  }
  .validate_precision_fit_for_targets(
    precision_fit,
    targets = restricted_set,
    p = p,
    name = "precision_fit"
  )

  if (!coordinate_restriction && is.null(C_standardized)) {
    stop("Internal restriction standardization failure.", call. = FALSE)
  }
  score_state <- if (reuse_full_score) {
    .subset_score_state(full_debias_state, restricted_set)
  } else {
    .make_score_state(
      X = prep$X,
      y = prep$y,
      beta_full = full_solver$beta,
      restricted_set = restricted_set,
      C_standardized = C_standardized,
      precision = precision_fit,
      x_scale = prep$x_scale,
      coordinate_restriction = coordinate_restriction
    )
  }

  if (test == "auto") {
    test <- if (q <= fixed_q_max) "wald" else "max"
  }
  if (test == "wald" && q > fixed_q_max) {
    stop(
      "Wald calibration is limited to q <= fixed_q_max. Use test = 'max' ",
      "for a large restriction block.",
      call. = FALSE
    )
  }
  inference <- if (test == "wald") {
    .wald_inference(score_state, target = t, alpha = alpha)
  } else {
    .max_inference(
      score_state = score_state,
      target = t,
      alpha = alpha,
      bootstrap_B = bootstrap_B,
      bootstrap_seed = bootstrap_seed,
      bootstrap_multipliers = bootstrap_multipliers,
      bootstrap_block_size = bootstrap_block_size
    )
  }

  beta_full_regularized <- full_solver$beta / prep$x_scale
  beta_submodel_regularized <- restricted_solver$beta / prep$x_scale
  beta_full <- if (endpoint == "debiased_models") {
    beta_full_debiased
  } else {
    beta_full_regularized
  }
  beta_submodel <- if (endpoint == "debiased_models") {
    beta_submodel_debiased
  } else {
    beta_submodel_regularized
  }
  shrinkage <- if (test == "wald") {
    .make_wald_shrinkage(
      beta_full,
      beta_submodel,
      statistic = inference$statistic,
      q = q,
      reject = inference$reject,
      epsilon = epsilon
    )
  } else {
    .make_max_shrinkage(
      beta_full,
      beta_submodel,
      statistic = inference$statistic,
      null_energy = inference$null_energy,
      reject = inference$reject,
      epsilon = epsilon
    )
  }

  all_beta <- list(
    full = beta_full,
    submodel = beta_submodel,
    restricted = beta_submodel,
    full_regularized = beta_full_regularized,
    submodel_regularized = beta_submodel_regularized,
    full_debiased = beta_full_debiased,
    submodel_debiased = beta_submodel_debiased,
    preliminary_test = shrinkage$preliminary_test,
    stein = shrinkage$stein,
    positive_part = shrinkage$positive_part,
    protected_positive_part = shrinkage$protected_positive_part
  )
  intercept <- lapply(all_beta, function(beta) {
    if (is.null(beta)) return(NULL)
    prep$y_center - sum(prep$x_center * beta)
  })

  restriction <- list(
    type = if (coordinate_restriction) "coordinate" else "linear",
    M = restricted_set,
    C = if (coordinate_restriction) NULL else C,
    t = t,
    zero_block = zero_block,
    construction = "projection_of_full",
    geometry = "standardized_euclidean"
  )
  output <- list(
    call = call,
    dimensions = c(
      n = n,
      p = p,
      p1 = p - m,
      p2 = m,
      m = m,
      q = q
    ),
    beta = all_beta,
    intercept = intercept,
    endpoint = endpoint,
    test = test,
    reject = inference$reject,
    alpha = alpha,
    inference = inference,
    shrinkage = shrinkage,
    full_solver = full_solver,
    restricted_solver = restricted_solver,
    restriction_projection = restricted_solver,
    restriction = restriction,
    precision = precision_fit,
    endpoint_precision = endpoint_precision,
    score_representation = score_state$direction_representation,
    preprocessing = prep[c("x_center", "x_scale", "y_center")],
    standardized_restriction = standardized_restriction
  )
  class(output) <- "hd_shrinkage_fit"
  output
}
