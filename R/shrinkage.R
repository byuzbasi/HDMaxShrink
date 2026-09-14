.make_wald_shrinkage <- function(
    beta_full,
    beta_restricted,
    statistic,
    q,
    reject,
    epsilon) {
  difference <- beta_full - beta_restricted
  preliminary_test <- if (reject) beta_full else beta_restricted
  if (q < 3L) {
    return(list(
      preliminary_test = preliminary_test,
      stein = NULL,
      positive_part = NULL,
      stein_weight = NA_real_,
      positive_weight = NA_real_,
      calibration = NA_real_,
      message = "Classical Wald-Stein estimators require q >= 3."
    ))
  }
  calibration <- q - 2
  statistic_safe <- max(statistic, epsilon)
  weight <- 1 - calibration / statistic_safe
  list(
    preliminary_test = preliminary_test,
    stein = beta_restricted + weight * difference,
    positive_part = beta_restricted + max(0, weight) * difference,
    stein_weight = weight,
    positive_weight = max(0, weight),
    calibration = calibration,
    message = NULL
  )
}

.make_max_shrinkage <- function(
    beta_full,
    beta_restricted,
    statistic,
    null_energy = NULL,
    reject,
    epsilon,
    shrinkage_calibration = NULL,
    shrinkage_calibration_type = NULL) {
  if (is.null(shrinkage_calibration)) {
    shrinkage_calibration <- null_energy
    shrinkage_calibration_type <-
      shrinkage_calibration_type %||% "second_moment"
  }
  if (length(shrinkage_calibration) != 1L ||
      !is.finite(shrinkage_calibration) || shrinkage_calibration <= 0) {
    stop(
      "The max-Stein calibration must be finite and positive.",
      call. = FALSE
    )
  }
  if (length(shrinkage_calibration_type) != 1L ||
      !is.character(shrinkage_calibration_type) ||
      is.na(shrinkage_calibration_type) ||
      !nzchar(shrinkage_calibration_type)) {
    stop("The max-Stein calibration type must be explicit.", call. = FALSE)
  }
  difference <- beta_full - beta_restricted
  statistic_squared <- max(statistic^2, epsilon)
  weight <- 1 - shrinkage_calibration / statistic_squared
  positive_weight <- max(0, weight)
  protected_weight <- if (reject) 1 else positive_weight
  list(
    type = "max_stein",
    preliminary_test = if (reject) beta_full else beta_restricted,
    stein = beta_restricted + weight * difference,
    positive_part = beta_restricted + positive_weight * difference,
    protected_positive_part =
      beta_restricted + protected_weight * difference,
    stein_weight = weight,
    positive_weight = positive_weight,
    protected_weight = protected_weight,
    calibration = shrinkage_calibration,
    calibration_type = shrinkage_calibration_type,
    shrinkage_calibration = shrinkage_calibration,
    shrinkage_calibration_type = shrinkage_calibration_type,
    message = paste0(
      "Max-calibrated Stein shrinkage is experimental; classical ",
      "James-Stein dominance is not asserted. The protected positive-part ",
      "estimator is reported separately and does not replace PS."
    )
  )
}
