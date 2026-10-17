# CPSS–SM shrinkage toward an all-predictor dual-Ridge endpoint

`fit_cpss_ridge_shrinkage()` combines an all-`X` Ridge full-model
endpoint with a fresh exact-null SVD submodel endpoint. A
conditional-Gaussian maximum partial-t test is calculated from the
exact-null residual, not from the biased Ridge estimator. The function
returns preliminary-test, Stein-type, positive-part, and
pretest-protected positive-part estimators.

## Usage

``` r
fit_cpss_ridge_shrinkage(
  X,
  y,
  core_set = NULL,
  selection = NULL,
  ridge_lambda = NULL,
  ridge_lambda_grid = NULL,
  selector_label = NULL,
  alpha = 0.05,
  gaussian_draws = NULL,
  bootstrap_B = 999L,
  bootstrap_seed = NULL,
  bootstrap_block_size = 256L,
  rank_tolerance = 1e-10,
  shrinkage_calibration = c("inverse_moment", "second_moment"),
  epsilon = 1e-10
)
```

## Arguments

- X:

  Numeric analysis-sample design matrix.

- y:

  Numeric analysis-sample response.

- core_set:

  One-based selected core indices, or a `cpss_core_selection` object.

- selection:

  Alternatively, an object returned by
  [`cpss_select_core()`](https://byuzbasi.github.io/HDMaxShrink/reference/cpss_select_core.md).
  Named features are aligned to `X` by identity.

- ridge_lambda:

  Optional fixed positive Ridge penalty. If omitted, GCV selects from
  `ridge_lambda_grid` using the existing all-`X` dual path.

- ridge_lambda_grid:

  Optional positive GCV grid.

- selector_label:

  Label describing the independently selected core.

- alpha:

  Conditional Monte Carlo test level.

- gaussian_draws:

  Optional reusable `n` by `bootstrap_B` matrix of iid standard Gaussian
  draws.

- bootstrap_B:

  Number of Gaussian draws when `gaussian_draws` is absent.

- bootstrap_seed:

  Optional reproducible Gaussian-draw seed.

- bootstrap_block_size:

  C++ calibration block size.

- rank_tolerance:

  Relative economy-SVD rank tolerance.

- shrinkage_calibration:

  Null calibration for the max-Stein weight.

- epsilon:

  Positive numerical floor for the squared statistic.

## Value

An `hd_shrinkage_fit` object with `FM`, `SM`, `PT`, `S`, `PS`, and `PPS`
endpoints.

## References

Shah, R. D. and Samworth, R. J. (2013). Variable selection with error
control: another look at stability selection. *Journal of the Royal
Statistical Society: Series B*, 75, 55–80.
[doi:10.1111/j.1467-9868.2011.01034.x](https://doi.org/10.1111/j.1467-9868.2011.01034.x)

Dufour, J.-M. (2006). Monte Carlo tests with nuisance parameters: a
general approach to finite-sample inference and nonstandard asymptotics.
*Journal of Econometrics*, 133, 443–477.
[doi:10.1016/j.jeconom.2005.06.007](https://doi.org/10.1016/j.jeconom.2005.06.007)

Yuzbasi, B., Arashi, M. and Ahmed, S. E. (2020). Shrinkage estimation
strategies in high-dimensional linear models. *International Statistical
Review*. [doi:10.1111/insr.12351](https://doi.org/10.1111/insr.12351)
