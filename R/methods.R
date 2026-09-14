.normalize_method <- function(method) {
  choices <- c(
    "full", "submodel", "restricted", "FM", "SM",
    "RIDGE", "MCP", "FSL", "SMSL", "FDSL", "SMDSL", "PT", "S", "PS",
    "PPS"
  )
  method <- match.arg(method, choices)
  switch(
    method,
    full = "full",
    submodel = "submodel",
    restricted = "restricted",
    FM = "full",
    SM = "submodel",
    RIDGE = "ridge",
    MCP = "mcp",
    FSL = "full_regularized",
    SMSL = "submodel_regularized",
    FDSL = "full_debiased",
    SMDSL = "submodel_debiased",
    PT = "preliminary_test",
    S = "stein",
    PS = "positive_part",
    PPS = "protected_positive_part"
  )
}

#' @export
print.hd_shrinkage_fit <- function(x, ...) {
  dimensions <- x$dimensions
  full_endpoint <- x$full_endpoint %||% switch(
    x$full_solver$construction,
    profiled_partial_square_root_lasso = "profiled_square_root_lasso",
    all_x_square_root_lasso = "all_x_square_root_lasso",
    "thresholded_partial_debiased"
  )
  heading <- if (identical(x$endpoint, "ridge_mcp")) {
    "Tmax-guided all-X Ridge--MCP shrinkage"
  } else if (identical(
    full_endpoint, "profiled_square_root_lasso"
  )) {
    "Tmax-guided profiled partial square-root-LASSO exact-null shrinkage"
  } else if (identical(full_endpoint, "all_x_square_root_lasso")) {
    "Tmax-guided all-X square-root-LASSO exact-null shrinkage"
  } else if (identical(x$endpoint, "partial_sqrt_lasso_null")) {
    "Tmax-guided thresholded partial-debiased exact-null shrinkage"
  } else {
    "Restriction-adaptive high-dimensional shrinkage"
  }
  cat(heading, "\n")
  cat(sprintf(
    paste0(
      "n=%d, p=%d, p1=%d, p2=%d, q=%d; endpoint=%s; ",
      "test=%s; reject=%s\n"
    ),
    dimensions["n"], dimensions["p"], dimensions["p1"],
    dimensions["p2"], dimensions["q"], x$endpoint, x$test, x$reject
  ))
  statistic_name <- if (x$test == "wald") {
    "W"
  } else if (identical(x$test, "max_partial_t")) {
    "max partial-t"
  } else {
    "Tmax"
  }
  cat(sprintf(
    "%s=%.5f; critical=%.5f; p-value=%.5f\n",
    statistic_name,
    x$inference$statistic,
    x$inference$critical_value,
    x$inference$p_value
  ))
  if (identical(x$endpoint, "ridge_mcp")) {
    cat(sprintf(
      paste0(
        "Ridge selection=%s; MCP selection=%s; MCP df=%d; ",
        "MCP convergence=%s\n"
      ),
      x$full_solver$selection,
      x$restricted_solver$selection,
      x$restricted_solver$degrees_freedom,
      x$restricted_solver$converged
    ))
  } else if (identical(x$endpoint, "partial_sqrt_lasso_null") &&
             identical(full_endpoint, "all_x_square_root_lasso")) {
    cat(sprintf(
      paste0(
        "full convergence=%s; FM=one-shot all-X square-root LASSO; ",
        "active all-X=%d; standardized intercept=0; ",
        "submodel=exact null SVD refit on X_A; core rank=%d; ",
        "restriction error=%.3e\n"
      ),
      x$full_solver$converged,
      length(x$full_solver$support_all),
      x$restricted_solver$rank,
      x$restricted_solver$restriction_violation
    ))
  } else if (identical(x$endpoint, "partial_sqrt_lasso_null") &&
             identical(full_endpoint, "profiled_square_root_lasso")) {
    cat(sprintf(
      paste0(
        "pilot convergence=%s; FM=profiled partial square-root LASSO; ",
        "active tested=%d; submodel=exact null SVD refit; core rank=%d; ",
        "restriction error=%.3e\n"
      ),
      x$full_solver$converged,
      length(x$full_solver$support_tested),
      x$restricted_solver$rank,
      x$restricted_solver$restriction_violation
    ))
  } else if (identical(x$endpoint, "partial_sqrt_lasso_null")) {
    cat(sprintf(
      paste0(
        "pilot convergence=%s; FM threshold eta=%.3g; retained=%d; ",
        "submodel=exact null SVD refit; core rank=%d; ",
        "restriction error=%.3e\n"
      ),
      x$full_solver$converged,
      x$threshold$eta,
      length(x$threshold$support_tested),
      x$restricted_solver$rank,
      x$restricted_solver$restriction_violation
    ))
  } else {
    cat(sprintf(
      paste0(
        "full convergence=%s; submodel=%s; ",
        "restriction error=%.3e\n"
      ),
      x$full_solver$converged,
      x$restriction$construction,
      x$restricted_solver$restriction_violation
    ))
  }
  if (x$test %in% c("max", "max_partial_t")) {
    calibration <- x$shrinkage$shrinkage_calibration %||%
      x$shrinkage$calibration
    calibration_type <- x$shrinkage$shrinkage_calibration_type %||%
      x$shrinkage$calibration_type %||% "second_moment"
    cat(sprintf(
      paste0(
        "max-Stein calibration (%s)=%.5f; PS weight=%.5f; ",
        "PPS weight=%.5f\n"
      ),
      calibration_type,
      calibration,
      x$shrinkage$positive_weight,
      x$shrinkage$protected_weight %||% NA_real_
    ))
  }
  if (!is.null(x$shrinkage$message)) {
    cat(x$shrinkage$message, "\n")
  }
  invisible(x)
}

#' Extract coefficients from a high-dimensional shrinkage fit
#'
#' @param object A fitted `hd_shrinkage_fit` object.
#' @param method Selected endpoint (`"FM"` or `"SM"`), adaptive estimator
#'   (`"PT"`, `"S"`, `"PS"`, or `"PPS"`), raw square-root lasso endpoint
#'   (`"FSL"` or
#'   `"SMSL"`), explicit debiased endpoint (`"FDSL"` or `"SMDSL"`) when
#'   available, or the
#'   all-`X` endpoints `"RIDGE"` and `"MCP"` from
#'   [fit_ridge_mcp_shrinkage()]. The legacy names `"full"` and `"restricted"`
#'   remain available.
#' @param ... Unused.
#'
#' @return A named coefficient vector including the intercept.
#' @export
coef.hd_shrinkage_fit <- function(
    object,
    method = c(
      "PS", "PPS", "FM", "SM", "PT", "S", "RIDGE", "MCP",
      "FSL", "SMSL", "FDSL", "SMDSL",
      "full", "submodel", "restricted"
    ),
    ...) {
  method <- match.arg(method)
  key <- .normalize_method(method)
  beta <- object$beta[[key]]
  intercept <- object$intercept[[key]]
  if (is.null(beta)) {
    stop("The requested estimator is undefined for this fit.", call. = FALSE)
  }
  result <- c(`(Intercept)` = intercept, beta)
  names(result)[-1L] <- paste0("x", seq_along(beta))
  result
}

#' Predict from a high-dimensional shrinkage fit
#'
#' @param object A fitted `hd_shrinkage_fit` object.
#' @param newdata Numeric matrix with the original design columns.
#' @param method Selected endpoint (`"FM"` or `"SM"`), adaptive estimator
#'   (`"PT"`, `"S"`, `"PS"`, or `"PPS"`), raw square-root lasso endpoint
#'   (`"FSL"` or
#'   `"SMSL"`), explicit debiased endpoint (`"FDSL"` or `"SMDSL"`) when
#'   available, or the
#'   all-`X` endpoints `"RIDGE"` and `"MCP"` from
#'   [fit_ridge_mcp_shrinkage()]. The legacy names `"full"` and `"restricted"`
#'   remain available.
#' @param ... Unused.
#'
#' @return Numeric predictions.
#' @export
predict.hd_shrinkage_fit <- function(
    object,
    newdata,
    method = c(
      "PS", "PPS", "FM", "SM", "PT", "S", "RIDGE", "MCP",
      "FSL", "SMSL", "FDSL", "SMDSL",
      "full", "submodel", "restricted"
    ),
    ...) {
  method <- match.arg(method)
  key <- .normalize_method(method)
  beta <- object$beta[[key]]
  intercept <- object$intercept[[key]]
  if (is.null(beta)) {
    stop("The requested estimator is undefined for this fit.", call. = FALSE)
  }
  newdata <- .validate_numeric_matrix(newdata, "newdata")
  if (ncol(newdata) != length(beta)) {
    stop("newdata must have the original p design columns.", call. = FALSE)
  }
  as.vector(intercept + newdata %*% beta)
}
