#' HDMaxShrink: test-guided high-dimensional shrinkage
#'
#' The package supports an independently learned CPSS core, a fresh exact-null
#' SVD submodel on that core, and an unrestricted all-predictor dual-Ridge full
#' endpoint. It also retains the all-predictor square-root-LASSO framework.
#' Both center and RMS-scale `X` and `y` and fix the standardized intercept at
#' zero. A conditional-Gaussian
#' maximum partial-t
#' test calibrates preliminary-test, Stein-type, positive-part, and separately
#' reported pretest-protected positive-part estimators. Economy-SVD and
#' compiled RcppArmadillo kernels
#' avoid OLS Gram inverses, `q` by `q` covariance inverses, and `p` by `p`
#' inverses. Earlier Ridge--MCP and restriction-projection modes remain
#' available for sensitivity analysis and reproducibility.
#'
#' @keywords internal
#' @useDynLib HDMaxShrink, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
