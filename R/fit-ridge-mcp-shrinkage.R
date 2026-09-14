.default_ridge_grid <- function() {
  exp(seq(log(1e-6), log(1e3), length.out = 100L))
}

.default_mcp_grid <- function(X, y, nlambda, lambda_min_ratio) {
  lambda_max <- max(abs(as.vector(crossprod(X, y)))) / nrow(X)
  if (!is.finite(lambda_max) || lambda_max <= .Machine$double.eps) {
    stop("The MCP lambda maximum is numerically zero.", call. = FALSE)
  }
  exp(seq(
    log(lambda_max),
    log(lambda_max * lambda_min_ratio),
    length.out = nlambda
  ))
}

.select_ridge_endpoint <- function(X, y, lambda, lambda_grid) {
  fixed <- !is.null(lambda)
  grid <- if (fixed) {
    .check_positive_scalar(lambda, "ridge_lambda")
  } else {
    lambda_grid %||% .default_ridge_grid()
  }
  grid <- sort(unique(as.numeric(grid)))
  if (!length(grid) || any(!is.finite(grid)) || any(grid <= 0)) {
    stop("ridge_lambda_grid must contain finite positive values.", call. = FALSE)
  }
  path <- cpp_ridge_dual_path(X, y, grid)
  selected <- if (fixed) 1L else which.min(path$gcv)
  list(
    beta = as.numeric(path$beta[, selected]),
    lambda = path$lambda[selected],
    selected_index = selected,
    selection = if (fixed) "fixed" else "GCV",
    rss = path$rss[selected],
    degrees_freedom = path$degrees_freedom[selected],
    criterion = path$gcv[selected],
    path = path[c(
      "lambda", "rss", "degrees_freedom", "gcv",
      "kernel_eigenvalues", "converged"
    )],
    converged = isTRUE(path$converged),
    iterations = 1L
  )
}

.mcp_ebic <- function(rss, degrees_freedom, n, p, gamma) {
  rss <- pmax(as.numeric(rss), n * .Machine$double.eps)
  degrees_freedom <- as.integer(degrees_freedom)
  n * log(rss / n) + degrees_freedom * log(n) +
    2 * gamma * lchoose(p, pmin(degrees_freedom, p))
}

.select_mcp_endpoint <- function(
    X,
    y,
    lambda,
    lambda_grid,
    gamma,
    ebic_gamma,
    nlambda,
    lambda_min_ratio,
    beta_init,
    max_iter,
    tol,
    zero_tol) {
  fixed <- !is.null(lambda)
  grid <- if (fixed) {
    .check_positive_scalar(lambda, "mcp_lambda")
  } else {
    lambda_grid %||% .default_mcp_grid(
      X, y, nlambda = nlambda, lambda_min_ratio = lambda_min_ratio
    )
  }
  grid <- sort(unique(as.numeric(grid)), decreasing = TRUE)
  if (!length(grid) || any(!is.finite(grid)) || any(grid <= 0)) {
    stop("mcp_lambda_grid must contain finite positive values.", call. = FALSE)
  }
  beta_init <- beta_init %||% numeric(0)
  path <- cpp_mcp_path(
    X = X,
    y = y,
    lambda = grid,
    gamma = gamma,
    beta_init = as.numeric(beta_init),
    max_iter = as.integer(max_iter),
    tol = tol,
    zero_tol = zero_tol
  )
  ebic <- .mcp_ebic(
    path$rss, path$degrees_freedom,
    n = nrow(X), p = ncol(X), gamma = ebic_gamma
  )
  eligible <- which(
    as.logical(path$converged) & is.finite(ebic) &
      path$degrees_freedom < nrow(X)
  )
  if (!length(eligible)) {
    eligible <- which(as.logical(path$converged) & is.finite(ebic))
  }
  if (!length(eligible)) {
    stop("No converged MCP path point is available.", call. = FALSE)
  }
  selected <- if (fixed) {
    if (!isTRUE(path$converged[1L])) {
      stop("The fixed MCP fit did not converge.", call. = FALSE)
    }
    1L
  } else {
    eligible[which.min(ebic[eligible])]
  }
  selected_beta <- as.numeric(path$beta[, selected])
  list(
    beta = selected_beta,
    lambda = path$lambda[selected],
    gamma = gamma,
    selected_index = selected,
    selection = if (fixed) "fixed" else "EBIC",
    ebic_gamma = ebic_gamma,
    ebic = ebic[selected],
    rss = path$rss[selected],
    objective = path$objective[selected],
    degrees_freedom = path$degrees_freedom[selected],
    support = which(abs(selected_beta) > zero_tol),
    converged = isTRUE(path$converged[selected]),
    iterations = path$iterations[selected],
    maximum_change = path$maximum_change[selected],
    path = list(
      lambda = path$lambda,
      rss = path$rss,
      objective = path$objective,
      degrees_freedom = path$degrees_freedom,
      converged = path$converged,
      iterations = path$iterations,
      maximum_change = path$maximum_change,
      ebic = ebic
    )
  )
}

#' Tmax-guided shrinkage between all-X Ridge and all-X MCP estimators
#'
#' `fit_ridge_mcp_shrinkage()` uses every design column in both endpoints. The
#' dense full-model endpoint is Ridge, computed through an `n` by `n` dual
#' system, and the sparse-model endpoint is an MCP fit whose support is selected
#' from the data. No true active set or prespecified zero block is imposed on
#' the MCP optimization. The coordinate set supplied in `tested_set` is used
#' only by an auxiliary debiased square-root-LASSO maximum test.
#'
#' @param X Numeric `n` by `p` design matrix.
#' @param y Numeric response vector.
#' @param tested_set One-based coordinates in the null hypothesis.
#' @param t Coordinatewise null target; defaults to zero.
#' @param ridge_lambda Optional fixed Ridge penalty.
#' @param ridge_lambda_grid Positive Ridge grid selected by GCV when
#'   `ridge_lambda` is `NULL`. GCV counts the unpenalized intercept in its
#'   effective degrees of freedom.
#' @param mcp_lambda Optional fixed MCP penalty.
#' @param mcp_lambda_grid Decreasing MCP grid selected by EBIC when
#'   `mcp_lambda` is `NULL`.
#' @param mcp_gamma MCP concavity parameter; the default is 3.
#' @param mcp_ebic_gamma EBIC combinatorial-penalty multiplier.
#' @param mcp_nlambda Number of automatically generated MCP penalties.
#' @param mcp_lambda_min_ratio Smallest automatic MCP penalty relative to the
#'   zero-solution penalty.
#' @param mcp_beta_init Optional original-scale MCP warm start. It never
#'   encodes or accepts an active set.
#' @param mcp_max_iter Maximum coordinate-descent sweeps per path point.
#' @param mcp_tol MCP convergence tolerance.
#' @param mcp_zero_tol Numerical support threshold.
#' @param test_lambda Square-root-LASSO penalty used only for the test pilot.
#' @param precision Precision-direction method for the test pilot.
#' @param lambda_node Nodewise-LASSO penalty when applicable.
#' @param theta_rows Optional user precision rows.
#' @param precision_fit Optional reusable precision fit for `tested_set`.
#' @param assume_independent Required for diagonal precision.
#' @param alpha Test level.
#' @param bootstrap_B Number of multiplier-bootstrap draws.
#' @param bootstrap_seed Optional bootstrap seed.
#' @param bootstrap_multipliers Optional reusable multiplier matrix.
#' @param bootstrap_block_size Bootstrap block size.
#' @param pilot_beta_init Optional original-scale square-root-LASSO warm start.
#' @param pilot_operator_norm Optional reusable design operator norm.
#' @param pilot_max_iter Maximum pilot iterations.
#' @param pilot_tol Pilot convergence tolerance.
#' @param pilot_check_every Pilot diagnostic interval.
#' @param epsilon Numerical floor in the Stein weight.
#'
#' @return An `hd_shrinkage_fit` object whose `FM` endpoint is Ridge and whose
#'   `SM` endpoint is the data-adaptive all-`X` MCP estimate. `PPS` is reported
#'   separately from `PS`: it equals `FM` after rejection and `PS` otherwise.
#'
#' @details
#' With centered, standardized data, the full endpoint is
#' \deqn{\hat\beta^{FM}=X^T(XX^T+n\lambda_R I_n)^{-1}y.}
#' The sparse endpoint minimizes
#' \deqn{\|y-Xb\|_2^2/(2n)+\sum_j p_{\lambda_M,\gamma}^{MCP}(|b_j|).}
#' The test pilot is deliberately separate from both risk endpoints so that
#' Ridge bias is not used to calibrate the high-dimensional null test.
#' Max-Stein risk dominance is not asserted.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(40 * 80), 40, 80)
#' beta <- c(1.2, -1, 0.8, rep(0, 77))
#' y <- drop(X %*% beta + rnorm(40))
#' fit <- fit_ridge_mcp_shrinkage(
#'   X, y, tested_set = 4:80,
#'   precision = "diagonal", assume_independent = TRUE,
#'   bootstrap_B = 99, mcp_nlambda = 20,
#'   pilot_max_iter = 5000, pilot_tol = 1e-4
#' )
#' coef(fit, method = "PS")
#' coef(fit, method = "PPS")
#'
#' @export
fit_ridge_mcp_shrinkage <- function(
    X,
    y,
    tested_set,
    t = NULL,
    ridge_lambda = NULL,
    ridge_lambda_grid = NULL,
    mcp_lambda = NULL,
    mcp_lambda_grid = NULL,
    mcp_gamma = 3,
    mcp_ebic_gamma = 0.5,
    mcp_nlambda = 60L,
    mcp_lambda_min_ratio = 0.01,
    mcp_beta_init = NULL,
    mcp_max_iter = 5000L,
    mcp_tol = 1e-7,
    mcp_zero_tol = 1e-8,
    test_lambda = NULL,
    precision = c("nodewise", "diagonal", "user"),
    lambda_node = NULL,
    theta_rows = NULL,
    precision_fit = NULL,
    assume_independent = FALSE,
    alpha = 0.05,
    bootstrap_B = 999L,
    bootstrap_seed = NULL,
    bootstrap_multipliers = NULL,
    bootstrap_block_size = 256L,
    pilot_beta_init = NULL,
    pilot_operator_norm = NULL,
    pilot_max_iter = 30000L,
    pilot_tol = 5e-6,
    pilot_check_every = 50L,
    epsilon = 1e-10) {
  call <- match.call()
  precision <- match.arg(precision)
  prep <- .prepare_design(X, y)
  n <- nrow(prep$X)
  p <- ncol(prep$X)
  tested_set <- .validate_indices(tested_set, p, "tested_set")
  q <- length(tested_set)
  t <- as.numeric(t %||% numeric(q))
  if (length(t) != q || any(!is.finite(t))) {
    stop("t must contain length(tested_set) finite entries.", call. = FALSE)
  }
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie strictly between zero and one.", call. = FALSE)
  }
  mcp_gamma <- .check_positive_scalar(mcp_gamma, "mcp_gamma")
  if (mcp_gamma <= 1) stop("mcp_gamma must exceed one.", call. = FALSE)
  mcp_ebic_gamma <- .check_positive_scalar(
    mcp_ebic_gamma, "mcp_ebic_gamma", allow_zero = TRUE
  )
  mcp_nlambda <- as.integer(mcp_nlambda)
  if (length(mcp_nlambda) != 1L || is.na(mcp_nlambda) || mcp_nlambda < 2L) {
    stop("mcp_nlambda must be an integer of at least two.", call. = FALSE)
  }
  mcp_lambda_min_ratio <- .check_positive_scalar(
    mcp_lambda_min_ratio, "mcp_lambda_min_ratio"
  )
  if (mcp_lambda_min_ratio >= 1) {
    stop("mcp_lambda_min_ratio must be below one.", call. = FALSE)
  }
  mcp_tol <- .check_positive_scalar(mcp_tol, "mcp_tol")
  mcp_zero_tol <- .check_positive_scalar(
    mcp_zero_tol, "mcp_zero_tol", allow_zero = TRUE
  )
  epsilon <- .check_positive_scalar(epsilon, "epsilon")

  ridge <- .select_ridge_endpoint(
    prep$X, prep$y,
    lambda = ridge_lambda,
    lambda_grid = ridge_lambda_grid
  )
  mcp_init_standardized <- if (is.null(mcp_beta_init)) {
    NULL
  } else {
    mcp_beta_init <- as.numeric(mcp_beta_init)
    if (length(mcp_beta_init) != p || any(!is.finite(mcp_beta_init))) {
      stop("mcp_beta_init must contain p finite values.", call. = FALSE)
    }
    mcp_beta_init * prep$x_scale
  }
  mcp <- .select_mcp_endpoint(
    prep$X, prep$y,
    lambda = mcp_lambda,
    lambda_grid = mcp_lambda_grid,
    gamma = mcp_gamma,
    ebic_gamma = mcp_ebic_gamma,
    nlambda = mcp_nlambda,
    lambda_min_ratio = mcp_lambda_min_ratio,
    beta_init = mcp_init_standardized,
    max_iter = mcp_max_iter,
    tol = mcp_tol,
    zero_tol = mcp_zero_tol
  )

  test_lambda <- test_lambda %||%
    (1.1 * sqrt(2 * log(2 * p) / n))
  .check_positive_scalar(test_lambda, "test_lambda")
  pilot_init_standardized <- if (is.null(pilot_beta_init)) {
    numeric(p)
  } else {
    pilot_beta_init <- as.numeric(pilot_beta_init)
    if (length(pilot_beta_init) != p || any(!is.finite(pilot_beta_init))) {
      stop("pilot_beta_init must contain p finite values.", call. = FALSE)
    }
    pilot_beta_init * prep$x_scale
  }
  pilot <- .sqrt_lasso_standardized(
    X = prep$X,
    y = prep$y,
    lambda = test_lambda,
    beta_init = pilot_init_standardized,
    max_iter = pilot_max_iter,
    tol = pilot_tol,
    check_every = pilot_check_every,
    operator_norm_value = pilot_operator_norm
  )
  if (!isTRUE(pilot$converged)) {
    warning("The square-root-LASSO test pilot reached its iteration limit.")
  }
  if (is.null(precision_fit)) {
    precision_fit <- estimate_precision_rows(
      X = prep$X,
      targets = tested_set,
      method = precision,
      lambda_node = lambda_node,
      theta_rows = theta_rows,
      assume_independent = assume_independent
    )
  }
  .validate_precision_fit_for_targets(
    precision_fit, targets = tested_set, p = p,
    name = "precision_fit"
  )
  score_state <- .make_score_state(
    X = prep$X,
    y = prep$y,
    beta_full = pilot$beta,
    restricted_set = tested_set,
    C_standardized = NULL,
    precision = precision_fit,
    x_scale = prep$x_scale,
    coordinate_restriction = TRUE
  )
  inference <- .max_inference(
    score_state = score_state,
    target = t,
    alpha = alpha,
    bootstrap_B = bootstrap_B,
    bootstrap_seed = bootstrap_seed,
    bootstrap_multipliers = bootstrap_multipliers,
    bootstrap_block_size = bootstrap_block_size
  )

  beta_ridge <- ridge$beta / prep$x_scale
  beta_mcp <- mcp$beta / prep$x_scale
  shrinkage <- .make_max_shrinkage(
    beta_full = beta_ridge,
    beta_restricted = beta_mcp,
    statistic = inference$statistic,
    null_energy = inference$null_energy,
    reject = inference$reject,
    epsilon = epsilon
  )
  beta <- list(
    full = beta_ridge,
    submodel = beta_mcp,
    restricted = beta_mcp,
    ridge = beta_ridge,
    mcp = beta_mcp,
    sparse = beta_mcp,
    full_regularized = NULL,
    submodel_regularized = NULL,
    full_debiased = NULL,
    submodel_debiased = NULL,
    test_pilot = pilot$beta / prep$x_scale,
    preliminary_test = shrinkage$preliminary_test,
    stein = shrinkage$stein,
    positive_part = shrinkage$positive_part,
    protected_positive_part = shrinkage$protected_positive_part
  )
  intercept <- lapply(beta, function(coefficient) {
    if (is.null(coefficient)) return(NULL)
    prep$y_center - sum(prep$x_center * coefficient)
  })

  ridge$beta <- NULL
  mcp$beta <- NULL
  mcp$support_original_scale <- mcp$support
  restriction_error <- max(abs(beta_mcp[tested_set] - t))
  output <- list(
    call = call,
    dimensions = c(
      n = n, p = p, p1 = p - q, p2 = q, m = q, q = q
    ),
    beta = beta,
    intercept = intercept,
    endpoint = "ridge_mcp",
    test = "max",
    reject = inference$reject,
    alpha = alpha,
    inference = inference,
    shrinkage = shrinkage,
    full_solver = ridge,
    restricted_solver = c(
      mcp,
      list(
        construction = "data_adaptive_MCP_all_X",
        restriction_violation = restriction_error
      )
    ),
    restriction = list(
      type = "coordinate_test_only",
      M = tested_set,
      C = NULL,
      t = t,
      zero_block = all(t == 0),
      construction = "data_adaptive_MCP_all_X",
      exact_endpoint_restriction = FALSE
    ),
    precision = precision_fit,
    test_pilot_solver = pilot,
    score_representation = score_state$direction_representation,
    preprocessing = prep[c("x_center", "x_scale", "y_center")]
  )
  class(output) <- "hd_shrinkage_fit"
  output
}
