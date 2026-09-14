.make_score_state <- function(
    X,
    y,
    beta_full,
    restricted_set,
    C_standardized,
    precision,
    x_scale,
    coordinate_restriction) {
  if (precision$method == "diagonal" && coordinate_restriction) {
    target_point <- beta_full[restricted_set] / x_scale[restricted_set]
    direction_scale <-
      precision$direction_scale / x_scale[restricted_set]
    state <- cpp_debiased_scores_diagonal(
      X = X,
      y = y,
      beta_full = beta_full,
      targets = restricted_set,
      direction_scale = direction_scale,
      target_point = target_point
    )
    state$direction_representation <- "diagonal"
    return(state)
  }

  theta_rows <- precision$theta_rows
  if (is.null(theta_rows)) {
    p <- ncol(X)
    theta_rows <- matrix(0, nrow = length(restricted_set), ncol = p)
    theta_rows[cbind(seq_along(restricted_set), restricted_set)] <-
      precision$direction_scale
  }
  if (coordinate_restriction) {
    directions <- sweep(
      theta_rows, 1L, x_scale[restricted_set], "/"
    )
    target_point <- beta_full[restricted_set] / x_scale[restricted_set]
  } else {
    directions <- C_standardized %*% theta_rows
    target_point <- as.vector(C_standardized %*% beta_full[restricted_set])
  }
  state <- cpp_debiased_scores_dense(
    X = X,
    y = y,
    beta_full = beta_full,
    target_point = target_point,
    directions = directions
  )
  state$direction_representation <- "dense"
  state$directions <- directions
  state
}

.validate_precision_fit_for_targets <- function(
    precision_fit,
    targets,
    p,
    name = "precision_fit") {
  if (!is.list(precision_fit) || is.null(precision_fit$method) ||
      is.null(precision_fit$targets) ||
      !identical(as.integer(precision_fit$targets), as.integer(targets))) {
    stop(
      name, " must be a compatible estimate_precision_rows() result.",
      call. = FALSE
    )
  }
  if (precision_fit$method == "diagonal") {
    if (length(precision_fit$direction_scale) != length(targets) ||
        any(!is.finite(precision_fit$direction_scale)) ||
        any(precision_fit$direction_scale <= 0)) {
      stop(name, " has invalid diagonal directions.", call. = FALSE)
    }
  } else {
    theta_rows <- precision_fit$theta_rows
    if (!is.matrix(theta_rows) || nrow(theta_rows) != length(targets) ||
        ncol(theta_rows) != p || any(!is.finite(theta_rows))) {
      stop(name, " has incompatible precision rows.", call. = FALSE)
    }
  }
  invisible(precision_fit)
}

.subset_precision_fit <- function(precision_fit, targets) {
  locations <- match(targets, precision_fit$targets)
  if (anyNA(locations)) {
    stop("The requested precision rows are unavailable.", call. = FALSE)
  }
  output <- precision_fit
  output$targets <- as.integer(targets)
  if (precision_fit$method == "diagonal") {
    output$direction_scale <- precision_fit$direction_scale[locations]
  } else {
    output$theta_rows <- precision_fit$theta_rows[locations, , drop = FALSE]
    if (!is.null(precision_fit$tau2)) {
      output$tau2 <- precision_fit$tau2[locations]
    }
  }
  output
}

.debiased_coordinate_model <- function(
    X,
    y,
    beta_standardized,
    x_scale,
    precision_fit) {
  p <- ncol(X)
  targets <- seq_len(p)
  .validate_precision_fit_for_targets(
    precision_fit, targets = targets, p = p,
    name = "model debiasing precision fit"
  )
  target_point <- beta_standardized / x_scale
  if (precision_fit$method == "diagonal") {
    state <- cpp_debiased_scores_diagonal(
      X = X,
      y = y,
      beta_full = beta_standardized,
      targets = targets,
      direction_scale = precision_fit$direction_scale / x_scale,
      target_point = target_point
    )
    state$direction_representation <- "diagonal"
    return(state)
  }
  directions <- sweep(precision_fit$theta_rows, 1L, x_scale, "/")
  state <- cpp_debiased_scores_dense(
    X = X,
    y = y,
    beta_full = beta_standardized,
    target_point = target_point,
    directions = directions
  )
  state$direction_representation <- "dense"
  state$directions <- directions
  state
}

.subset_score_state <- function(score_state, targets) {
  output <- list(
    theta_tilde = score_state$theta_tilde[targets],
    residual = score_state$residual,
    psi_centered = score_state$psi_centered[, targets, drop = FALSE],
    variance = score_state$variance[targets],
    direction_representation = score_state$direction_representation
  )
  if (!is.null(score_state$directions)) {
    output$directions <- score_state$directions[targets, , drop = FALSE]
  }
  output
}

.wald_inference <- function(score_state, target, alpha) {
  psi <- score_state$psi_centered
  n <- nrow(psi)
  q <- ncol(psi)
  covariance <- crossprod(psi) / n
  covariance <- (covariance + t(covariance)) / 2
  if (any(!is.finite(covariance)) || any(diag(covariance) <= 0)) {
    stop(
      "The small-q contrast covariance is not positive on its diagonal.",
      call. = FALSE
    )
  }
  chol_covariance <- tryCatch(chol(covariance), error = function(e) NULL)
  if (is.null(chol_covariance)) {
    stop(
      "The small-q contrast covariance is singular; no ridge or diagonal ",
      "loading was applied.",
      call. = FALSE
    )
  }
  difference <- score_state$theta_tilde - target
  whitened <- backsolve(chol_covariance, difference, transpose = TRUE)
  statistic <- n * sum(whitened^2)
  critical <- stats::qchisq(1 - alpha, df = q)
  list(
    type = "wald",
    statistic = statistic,
    W = statistic,
    df = q,
    covariance = covariance,
    standard_error = sqrt(diag(covariance)),
    critical_value = critical,
    p_value = stats::pchisq(statistic, df = q, lower.tail = FALSE),
    reject = statistic > critical,
    theta_tilde = score_state$theta_tilde,
    target = target,
    psi_centered = psi
  )
}

#' Multiplier-bootstrap maximum test from debiased influence scores
#'
#' The function requires only coordinatewise score variances. It does not form
#' or invert a `q` by `q` covariance matrix and therefore permits `q > n`.
#'
#' @param theta_tilde Debiased target estimates.
#' @param target Null target vector.
#' @param psi_centered Centered `n` by `q` influence-score matrix.
#' @param alpha Test level.
#' @param bootstrap_B Number of Gaussian multiplier draws.
#' @param bootstrap_seed Optional reproducible seed. The caller's RNG state is
#'   restored after multiplier generation.
#' @param multipliers Optional fixed `n` by `bootstrap_B` multiplier matrix.
#' @param block_size Number of draws processed together by the C++ kernel.
#'
#' @return A list containing the observed maximum, bootstrap distribution,
#'   critical value, p-value, separate null-moment calibrations, and the
#'   backward-compatible `null_energy` alias for the null second moment.
#' @export
max_test_hd <- function(
    theta_tilde,
    target,
    psi_centered,
    alpha = 0.05,
    bootstrap_B = 999L,
    bootstrap_seed = NULL,
    multipliers = NULL,
    block_size = 256L) {
  theta_tilde <- as.numeric(theta_tilde)
  target <- as.numeric(target)
  psi_centered <- .validate_numeric_matrix(psi_centered, "psi_centered")
  n <- nrow(psi_centered)
  q <- ncol(psi_centered)
  if (length(theta_tilde) != q || length(target) != q ||
      any(!is.finite(theta_tilde)) || any(!is.finite(target))) {
    stop(
      "theta_tilde and target must contain q finite entries.",
      call. = FALSE
    )
  }
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie strictly between zero and one.", call. = FALSE)
  }

  variance <- colMeans(psi_centered^2)
  if (any(!is.finite(variance)) || any(variance <= 0)) {
    stop("Every score variance must be finite and positive.", call. = FALSE)
  }
  standard_error <- sqrt(variance)
  observed <- max(
    sqrt(n) * abs(theta_tilde - target) / standard_error
  )

  if (is.null(multipliers)) {
    bootstrap_B <- as.integer(bootstrap_B)
    if (length(bootstrap_B) != 1L || is.na(bootstrap_B) ||
        bootstrap_B < 99L) {
      stop("bootstrap_B must be an integer of at least 99.", call. = FALSE)
    }
    multipliers <- .with_seed(bootstrap_seed, function() {
      matrix(stats::rnorm(n * bootstrap_B), nrow = n, ncol = bootstrap_B)
    })
  } else {
    multipliers <- .validate_numeric_matrix(multipliers, "multipliers")
    if (nrow(multipliers) != n || ncol(multipliers) < 99L) {
      stop(
        "multipliers must have nrow(psi_centered) rows and at least 99 columns.",
        call. = FALSE
      )
    }
    bootstrap_B <- ncol(multipliers)
  }

  bootstrap_statistics <- cpp_multiplier_max(
    psi_centered = psi_centered,
    standard_error = standard_error,
    multipliers = multipliers,
    block_size = as.integer(block_size)
  )
  critical <- as.numeric(stats::quantile(
    bootstrap_statistics,
    probs = 1 - alpha,
    names = FALSE,
    type = 8
  ))
  null_moments <- .max_null_moments(bootstrap_statistics)

  list(
    type = "max",
    subtype = "residual_multiplier_debiased_max",
    statistic = observed,
    T_max = observed,
    score_statistic = observed,
    max_partial_t = NA_real_,
    critical_value = critical,
    p_value = (1 + sum(bootstrap_statistics >= observed)) /
      (bootstrap_B + 1),
    reject = observed > critical,
    standard_error = standard_error,
    variance = variance,
    bootstrap_statistics = bootstrap_statistics,
    null_second_moment = null_moments$second_moment,
    null_inverse_moment_calibration = null_moments$inverse_moment,
    shrinkage_calibration = null_moments$second_moment,
    shrinkage_calibration_type = "second_moment",
    null_energy = null_moments$second_moment,
    bootstrap_B = bootstrap_B,
    alpha = alpha,
    theta_tilde = theta_tilde,
    target = target
  )
}

.max_null_moments <- function(bootstrap_statistics) {
  bootstrap_statistics <- as.numeric(bootstrap_statistics)
  squared <- bootstrap_statistics^2
  if (!length(squared) || any(!is.finite(squared)) || any(squared <= 0)) {
    stop(
      "Bootstrap squared maxima must be finite and strictly positive.",
      call. = FALSE
    )
  }
  list(
    second_moment = mean(squared),
    inverse_moment = 1 / mean(1 / squared)
  )
}

.set_max_shrinkage_calibration <- function(
    inference,
    calibration_type = c("inverse_moment", "second_moment")) {
  calibration_type <- match.arg(calibration_type)
  moments <- .max_null_moments(inference$bootstrap_statistics)
  selected <- switch(
    calibration_type,
    inverse_moment = moments$inverse_moment,
    second_moment = moments$second_moment
  )
  inference$null_second_moment <- moments$second_moment
  inference$null_inverse_moment_calibration <- moments$inverse_moment
  inference$shrinkage_calibration <- selected
  inference$shrinkage_calibration_type <- calibration_type
  # Backward-compatible alias with its original, unambiguous meaning.
  inference$null_energy <- moments$second_moment
  inference
}

.monte_carlo_max_decision <- function(observed, bootstrap_statistics, alpha) {
  bootstrap_statistics <- as.numeric(bootstrap_statistics)
  draws <- length(bootstrap_statistics)
  if (draws < 99L || any(!is.finite(bootstrap_statistics)) ||
      length(observed) != 1L || !is.finite(observed)) {
    stop("Invalid Monte Carlo maximum-test inputs.", call. = FALSE)
  }
  p_value <- (1 + sum(bootstrap_statistics >= observed)) / (draws + 1)
  rejection_ranks <- floor(alpha * (draws + 1))
  if (rejection_ranks < 1L) {
    critical <- Inf
  } else {
    critical_index <- draws + 1L - rejection_ranks
    critical <- sort.int(
      bootstrap_statistics, partial = critical_index
    )[critical_index]
  }
  list(
    critical_value = critical,
    p_value = p_value,
    reject = p_value <= alpha
  )
}

.max_partial_t_transform <- function(score_statistic, residual_df) {
  score_statistic <- as.numeric(score_statistic)
  residual_df <- as.numeric(residual_df)
  if (length(residual_df) != 1L || !is.finite(residual_df) ||
      residual_df <= 1 || any(!is.finite(score_statistic)) ||
      any(score_statistic < 0)) {
    stop("Invalid max partial-t transformation inputs.", call. = FALSE)
  }
  statistic_squared <- score_statistic^2
  denominator <- residual_df - statistic_squared
  numerical_tolerance <- sqrt(.Machine$double.eps) * residual_df
  if (any(denominator < -numerical_tolerance)) {
    stop(
      "A score statistic exceeds its residual degrees-of-freedom bound.",
      call. = FALSE
    )
  }
  denominator <- pmax(denominator, .Machine$double.eps * residual_df)
  sqrt((residual_df - 1) * statistic_squared / denominator)
}

.as_max_partial_t_inference <- function(inference) {
  score_statistic <- as.numeric(inference$statistic)
  score_critical_value <- as.numeric(inference$critical_value)
  score_bootstrap_statistics <- as.numeric(inference$bootstrap_statistics)
  max_partial_t <- .max_partial_t_transform(
    score_statistic, inference$residual_df
  )
  bootstrap_partial_t <- .max_partial_t_transform(
    score_bootstrap_statistics, inference$residual_df
  )
  decision <- .monte_carlo_max_decision(
    observed = max_partial_t,
    bootstrap_statistics = bootstrap_partial_t,
    alpha = inference$alpha
  )
  if (!isTRUE(all.equal(decision$p_value, inference$p_value, tolerance = 0)) ||
      !identical(decision$reject, inference$reject)) {
    stop(
      "The monotone max partial-t transformation changed the test decision.",
      call. = FALSE
    )
  }
  inference$score_statistic <- score_statistic
  inference$score_critical_value <- score_critical_value
  inference$score_bootstrap_statistics <- score_bootstrap_statistics
  inference$max_partial_t <- max_partial_t
  inference$max_partial_t_critical_value <- decision$critical_value
  inference$max_partial_t_bootstrap_statistics <- bootstrap_partial_t
  inference$statistic <- max_partial_t
  inference$T_max <- max_partial_t
  inference$critical_value <- decision$critical_value
  inference$p_value <- decision$p_value
  inference$reject <- decision$reject
  inference$bootstrap_statistics <- bootstrap_partial_t
  inference$type <- "gaussian_max_partial_t"
  inference$subtype <- "conditional_gaussian_exact_null_max_partial_t"
  inference$statistic_scale <- "max_partial_t"
  inference$transformation <- paste0(
    "sqrt((df - 1) * T_score^2 / (df - T_score^2))"
  )
  inference
}

.gaussian_score_max_inference <- function(
    null_state,
    geometry,
    alpha = 0.05,
    bootstrap_B = 999L,
    bootstrap_seed = NULL,
    gaussian_draws = NULL,
    bootstrap_block_size = 256L) {
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie strictly between zero and one.", call. = FALSE)
  }
  n <- nrow(geometry$residualized_tested)
  if (is.null(gaussian_draws)) {
    bootstrap_B <- as.integer(bootstrap_B)
    if (length(bootstrap_B) != 1L || is.na(bootstrap_B) ||
        bootstrap_B < 99L) {
      stop("bootstrap_B must be an integer of at least 99.", call. = FALSE)
    }
    gaussian_draws <- .with_seed(bootstrap_seed, function() {
      matrix(stats::rnorm(n * bootstrap_B), nrow = n)
    })
  } else {
    gaussian_draws <- .validate_numeric_matrix(
      gaussian_draws, "gaussian_draws"
    )
    if (nrow(gaussian_draws) != n || ncol(gaussian_draws) < 99L) {
      stop(
        "gaussian_draws must have n rows and at least 99 columns.",
        call. = FALSE
      )
    }
    bootstrap_B <- ncol(gaussian_draws)
  }

  result <- cpp_gaussian_score_max(
    null_residual = null_state$residual,
    residualized_tested = geometry$residualized_tested,
    core_basis = geometry$orthonormal_basis,
    gaussian_draws = gaussian_draws,
    block_size = as.integer(bootstrap_block_size)
  )
  result$score <- as.numeric(result$score)
  result$direction_norm <- as.numeric(result$direction_norm)
  result$bootstrap_statistics <- as.numeric(result$bootstrap_statistics)
  decision <- .monte_carlo_max_decision(
    observed = result$statistic,
    bootstrap_statistics = result$bootstrap_statistics,
    alpha = alpha
  )
  result[names(decision)] <- decision
  result$type <- "gaussian_score_max"
  result$subtype <- "conditional_gaussian_exact_null_score"
  result$bootstrap_B <- bootstrap_B
  result$alpha <- alpha
  result$score_statistic <- result$statistic
  result$max_partial_t <- NA_real_
  result <- .set_max_shrinkage_calibration(
    result, calibration_type = "second_moment"
  )
  result$conditional_exactness <- paste0(
    "Finite-sample conditional Monte Carlo calibration under homoskedastic ",
    "Gaussian errors and a prespecified centered core."
  )
  result
}

.max_inference <- function(
    score_state,
    target,
    alpha,
    bootstrap_B,
    bootstrap_seed,
    bootstrap_multipliers,
    bootstrap_block_size) {
  result <- max_test_hd(
    theta_tilde = score_state$theta_tilde,
    target = target,
    psi_centered = score_state$psi_centered,
    alpha = alpha,
    bootstrap_B = bootstrap_B,
    bootstrap_seed = bootstrap_seed,
    multipliers = bootstrap_multipliers,
    block_size = bootstrap_block_size
  )
  result$psi_centered <- score_state$psi_centered
  result
}
