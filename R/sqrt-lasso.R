.sqrt_lasso_standardized <- function(
    X,
    y,
    lambda,
    A = NULL,
    target = NULL,
    penalty_factor = NULL,
    beta_init = NULL,
    max_iter = 30000L,
    tol = 5e-6,
    check_every = 50L,
    verbose = FALSE,
    operator_norm_value = NULL) {
  X <- .validate_numeric_matrix(X, "X")
  y <- .validate_response(y, nrow(X))
  p <- ncol(X)

  if (is.null(A)) {
    A <- matrix(0, nrow = 0L, ncol = p)
    target <- numeric(0)
  } else {
    A <- as.matrix(A)
    storage.mode(A) <- "double"
    if (ncol(A) != p || any(!is.finite(A))) {
      stop("A must be finite and have ncol(X) columns.", call. = FALSE)
    }
    target <- as.numeric(target)
    if (length(target) != nrow(A) || any(!is.finite(target))) {
      stop(
        "target must be finite and have nrow(A) entries.",
        call. = FALSE
      )
    }
  }
  penalty_factor <- penalty_factor %||% rep(1, p)
  beta_init <- beta_init %||% numeric(p)
  operator_norm_value <- operator_norm_value %||% NA_real_

  cpp_sqrt_lasso_pd(
    X = X,
    y = y,
    lambda = as.numeric(lambda),
    A = A,
    target = target,
    penalty_factor = as.numeric(penalty_factor),
    beta_init = as.numeric(beta_init),
    max_iter = as.integer(max_iter),
    tol = as.numeric(tol),
    check_every = as.integer(check_every),
    verbose = isTRUE(verbose),
    operator_norm_value = as.numeric(operator_norm_value)
  )
}

#' Fit a square-root lasso model with an RcppArmadillo solver
#'
#' This lower-level function fits an unconstrained square-root lasso after
#' centering and scaling the design. It never forms or inverts a Gram matrix.
#'
#' @param X Numeric design matrix.
#' @param y Numeric response vector.
#' @param lambda Positive penalty. The default is
#'   `1.1 * sqrt(2 * log(2 * p) / n)`.
#' @param penalty_factor Nonnegative penalty multiplier for every coefficient.
#' @param max_iter Maximum primal-dual iterations.
#' @param tol Convergence tolerance.
#' @param check_every Iteration interval for KKT checks.
#' @param verbose Print periodic solver diagnostics.
#'
#' @return A list containing coefficients, intercept, and solver diagnostics.
#' @export
sqrt_lasso_hd <- function(
    X,
    y,
    lambda = NULL,
    penalty_factor = NULL,
    max_iter = 30000L,
    tol = 5e-6,
    check_every = 50L,
    verbose = FALSE) {
  prep <- .prepare_design(X, y)
  n <- nrow(prep$X)
  p <- ncol(prep$X)
  lambda <- lambda %||% (1.1 * sqrt(2 * log(2 * p) / n))
  penalty_factor <- penalty_factor %||% rep(1, p)
  fit <- .sqrt_lasso_standardized(
    prep$X,
    prep$y,
    lambda = lambda,
    penalty_factor = penalty_factor,
    max_iter = max_iter,
    tol = tol,
    check_every = check_every,
    verbose = verbose
  )
  fit$beta_standardized <- fit$beta
  fit$beta <- fit$beta / prep$x_scale
  fit$intercept <- prep$y_center - sum(prep$x_center * fit$beta)
  fit$preprocessing <- prep[c("x_center", "x_scale", "y_center")]
  class(fit) <- c("hd_sqrt_lasso_fit", "list")
  fit
}
