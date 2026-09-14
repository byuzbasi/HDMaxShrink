`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.validate_numeric_matrix <- function(x, name) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (!length(x) || any(!is.finite(x))) {
    stop(name, " must be a non-empty finite numeric matrix.", call. = FALSE)
  }
  x
}

.validate_response <- function(y, n) {
  y <- as.numeric(y)
  if (length(y) != n || any(!is.finite(y))) {
    stop("y must be finite and have nrow(X) entries.", call. = FALSE)
  }
  y
}

.validate_indices <- function(index, p, name = "restricted_set") {
  index <- as.integer(index)
  if (!length(index) || anyNA(index) || any(index < 1L) ||
      any(index > p) || anyDuplicated(index)) {
    stop(
      name, " must contain unique one-based indices in 1,...,p.",
      call. = FALSE
    )
  }
  index
}

.check_positive_scalar <- function(x, name, allow_zero = FALSE) {
  x <- as.numeric(x)
  lower_ok <- if (allow_zero) x >= 0 else x > 0
  if (length(x) != 1L || !is.finite(x) || !lower_ok) {
    qualifier <- if (allow_zero) "nonnegative" else "positive"
    stop(name, " must be one finite ", qualifier, " number.", call. = FALSE)
  }
  x
}

.prepare_design <- function(X, y, scale_y = FALSE) {
  X <- .validate_numeric_matrix(X, "X")
  y <- .validate_response(y, nrow(X))
  scale_y <- isTRUE(scale_y)
  x_center <- colMeans(X)
  X_centered <- sweep(X, 2L, x_center, "-")
  x_scale <- sqrt(colMeans(X_centered^2))
  if (any(!is.finite(x_scale)) ||
      any(x_scale <= .Machine$double.eps^0.5)) {
    stop(
      "Every design column must have positive empirical variance.",
      call. = FALSE
    )
  }
  y_center <- mean(y)
  y_centered <- y - y_center
  y_scale <- if (scale_y) sqrt(mean(y_centered^2)) else 1
  if (!is.finite(y_scale) || y_scale <= .Machine$double.eps^0.5) {
    stop("The response must have positive empirical variance.", call. = FALSE)
  }
  list(
    X = sweep(X_centered, 2L, x_scale, "/"),
    y = y_centered / y_scale,
    x_center = x_center,
    x_scale = x_scale,
    y_center = y_center,
    y_scale = y_scale,
    standardized_intercept = 0
  )
}

.with_seed <- function(seed, callback) {
  if (is.null(seed)) {
    return(callback())
  }
  seed <- as.integer(seed)
  if (length(seed) != 1L || is.na(seed)) {
    stop("seed must be NULL or one finite integer.", call. = FALSE)
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  callback()
}

.is_identity_matrix <- function(C, tolerance = 1e-12) {
  nrow(C) == ncol(C) &&
    max(abs(C - diag(nrow(C)))) <= tolerance
}
