.project_full_standardized <- function(
    beta_full,
    restricted_set,
    coordinate_restriction,
    C_standardized,
    target) {
  beta_full <- as.numeric(beta_full)
  target <- as.numeric(target)

  if (coordinate_restriction) {
    beta_projected <- beta_full
    beta_projected[restricted_set] <- target
    return(list(
      beta = beta_projected,
      standardized_restriction = NULL,
      restriction_violation = 0
    ))
  }

  p <- length(beta_full)
  q <- nrow(C_standardized)
  A <- matrix(0, nrow = q, ncol = p)
  A[, restricted_set] <- C_standardized
  gram <- tcrossprod(C_standardized)
  chol_gram <- tryCatch(chol(gram), error = function(error) NULL)
  if (is.null(chol_gram)) {
    stop(
      "The standardized restriction geometry is singular; no ridge or ",
      "diagonal loading was applied.",
      call. = FALSE
    )
  }
  discrepancy <- as.vector(A %*% beta_full - target)
  multiplier <- backsolve(
    chol_gram,
    forwardsolve(t(chol_gram), discrepancy)
  )
  beta_projected <- beta_full - as.vector(crossprod(A, multiplier))
  violation <- max(abs(as.vector(A %*% beta_projected) - target))

  list(
    beta = beta_projected,
    standardized_restriction = A,
    restriction_violation = violation
  )
}

.make_projection_state <- function(
    full_solver,
    beta_projected,
    X,
    y,
    lambda,
    penalty_factor,
    restriction_violation) {
  residual <- as.vector(y - X %*% beta_projected)
  objective <- sqrt(sum(residual^2)) / sqrt(length(y)) +
    lambda * sum(penalty_factor * abs(beta_projected))

  list(
    beta = beta_projected,
    converged = isTRUE(full_solver$converged),
    iterations = 0L,
    objective = objective,
    relative_change = 0,
    kkt_residual = NA_real_,
    restriction_violation = restriction_violation,
    operator_norm = NA_real_,
    lambda = lambda,
    active_set = which(abs(beta_projected) > 0),
    construction = "projection_of_full",
    full_solver_converged = isTRUE(full_solver$converged)
  )
}
