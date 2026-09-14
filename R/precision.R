#' Estimate target rows of a high-dimensional precision matrix
#'
#' Nodewise lasso is the general option. The diagonal option is intentionally
#' restricted to designs whose independence/diagonal-covariance structure is
#' assumed in advance.
#'
#' @param X Centered and standardized numeric design matrix.
#' @param targets One-based target-column indices.
#' @param method One of `"nodewise"`, `"diagonal"`, or `"user"`.
#' @param lambda_node Nodewise lasso penalty.
#' @param theta_rows User-supplied rows when `method = "user"`.
#' @param assume_independent Must be `TRUE` for diagonal precision.
#' @param glmnet_thresh Convergence threshold passed to [glmnet::glmnet()].
#'
#' @return A precision-direction specification.
#' @export
estimate_precision_rows <- function(
    X,
    targets,
    method = c("nodewise", "diagonal", "user"),
    lambda_node = NULL,
    theta_rows = NULL,
    assume_independent = FALSE,
    glmnet_thresh = 1e-10) {
  method <- match.arg(method)
  X <- .validate_numeric_matrix(X, "X")
  n <- nrow(X)
  p <- ncol(X)
  targets <- .validate_indices(targets, p, "targets")
  q <- length(targets)

  if (method == "diagonal") {
    if (!isTRUE(assume_independent)) {
      stop(
        "precision = 'diagonal' requires assume_independent = TRUE.",
        call. = FALSE
      )
    }
    empirical_variance <- colMeans(X[, targets, drop = FALSE]^2)
    if (any(!is.finite(empirical_variance)) ||
        any(empirical_variance <= .Machine$double.eps^0.5)) {
      stop("Target variances must be finite and positive.", call. = FALSE)
    }
    return(list(
      method = method,
      targets = targets,
      direction_scale = 1 / empirical_variance,
      theta_rows = NULL,
      lambda = NA_real_
    ))
  }

  if (method == "user") {
    theta_rows <- .validate_numeric_matrix(theta_rows, "theta_rows")
    if (nrow(theta_rows) != q || ncol(theta_rows) != p) {
      stop(
        "theta_rows must have length(targets) rows and ncol(X) columns.",
        call. = FALSE
      )
    }
    return(list(
      method = method,
      targets = targets,
      direction_scale = NULL,
      theta_rows = theta_rows,
      lambda = NA_real_
    ))
  }

  lambda_node <- lambda_node %||%
    (1.1 * sqrt(2 * log(2 * p) / n))
  .check_positive_scalar(lambda_node, "lambda_node")
  theta_out <- matrix(0, nrow = q, ncol = p)
  tau2 <- numeric(q)
  all_columns <- seq_len(p)

  for (k in seq_along(targets)) {
    j <- targets[k]
    minus_j <- all_columns[all_columns != j]
    fit <- glmnet::glmnet(
      x = X[, minus_j, drop = FALSE],
      y = X[, j],
      family = "gaussian",
      alpha = 1,
      lambda = lambda_node,
      intercept = FALSE,
      standardize = FALSE,
      thresh = glmnet_thresh,
      maxit = 100000L
    )
    gamma <- as.numeric(stats::coef(fit, s = lambda_node))[-1L]
    residual <- X[, j] -
      as.vector(X[, minus_j, drop = FALSE] %*% gamma)
    tau2[k] <- mean(residual^2) + lambda_node * sum(abs(gamma))
    if (!is.finite(tau2[k]) ||
        tau2[k] <= .Machine$double.eps^0.5) {
      stop(
        "Nodewise residual scale is nonpositive for column ", j, ".",
        call. = FALSE
      )
    }
    direction <- numeric(p)
    direction[j] <- 1
    direction[minus_j] <- -gamma
    theta_out[k, ] <- direction / tau2[k]
  }

  list(
    method = method,
    targets = targets,
    direction_scale = NULL,
    theta_rows = theta_out,
    tau2 = tau2,
    lambda = lambda_node
  )
}
