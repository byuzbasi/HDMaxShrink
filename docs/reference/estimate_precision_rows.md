# Estimate target rows of a high-dimensional precision matrix

Nodewise lasso is the general option. The diagonal option is
intentionally restricted to designs whose
independence/diagonal-covariance structure is assumed in advance.

## Usage

``` r
estimate_precision_rows(
  X,
  targets,
  method = c("nodewise", "diagonal", "user"),
  lambda_node = NULL,
  theta_rows = NULL,
  assume_independent = FALSE,
  glmnet_thresh = 1e-10
)
```

## Arguments

- X:

  Centered and standardized numeric design matrix.

- targets:

  One-based target-column indices.

- method:

  One of `"nodewise"`, `"diagonal"`, or `"user"`.

- lambda_node:

  Nodewise lasso penalty.

- theta_rows:

  User-supplied rows when `method = "user"`.

- assume_independent:

  Must be `TRUE` for diagonal precision.

- glmnet_thresh:

  Convergence threshold passed to
  [`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html).

## Value

A precision-direction specification.
