# RcppArmadillo primal-dual square-root lasso solver

RcppArmadillo primal-dual square-root lasso solver

## Usage

``` r
cpp_sqrt_lasso_pd(
  X,
  y,
  lambda,
  A,
  target,
  penalty_factor,
  beta_init,
  max_iter = 30000L,
  tol = 5e-06,
  check_every = 50L,
  verbose = FALSE,
  operator_norm_value = NA_real_
)
```
