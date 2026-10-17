# Fit a square-root lasso model with an RcppArmadillo solver

This lower-level function fits an unconstrained square-root lasso after
centering and scaling the design. It never forms or inverts a Gram
matrix.

## Usage

``` r
sqrt_lasso_hd(
  X,
  y,
  lambda = NULL,
  penalty_factor = NULL,
  max_iter = 30000L,
  tol = 5e-06,
  check_every = 50L,
  verbose = FALSE
)
```

## Arguments

- X:

  Numeric design matrix.

- y:

  Numeric response vector.

- lambda:

  Positive penalty. The default is `1.1 * sqrt(2 * log(2 * p) / n)`.

- penalty_factor:

  Nonnegative penalty multiplier for every coefficient.

- max_iter:

  Maximum primal-dual iterations.

- tol:

  Convergence tolerance.

- check_every:

  Iteration interval for KKT checks.

- verbose:

  Print periodic solver diagnostics.

## Value

A list containing coefficients, intercept, and solver diagnostics.
