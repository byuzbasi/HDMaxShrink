# One-step scores for the profiled partial square-root LASSO

The tested design has already been residualized against the prespecified
core. A diagonal population precision is used, which is appropriate for
the Gaussian-iid design studied in the accompanying simulations. No
tested-block covariance matrix is formed or inverted.

## Usage

``` r
cpp_partial_debiased_scores(y_residualized, residualized_tested, beta_tested)
```
