# Tmax-guided shrinkage between all-X Ridge and all-X MCP estimators

`fit_ridge_mcp_shrinkage()` uses every design column in both endpoints.
The dense full-model endpoint is Ridge, computed through an `n` by `n`
dual system, and the sparse-model endpoint is an MCP fit whose support
is selected from the data. No true active set or prespecified zero block
is imposed on the MCP optimization. The coordinate set supplied in
`tested_set` is used only by an auxiliary debiased square-root-LASSO
maximum test.

## Usage

``` r
fit_ridge_mcp_shrinkage(
  X,
  y,
  tested_set,
  t = NULL,
  ridge_lambda = NULL,
  ridge_lambda_grid = NULL,
  mcp_lambda = NULL,
  mcp_lambda_grid = NULL,
  mcp_gamma = 3,
  mcp_ebic_gamma = 0.5,
  mcp_nlambda = 60L,
  mcp_lambda_min_ratio = 0.01,
  mcp_beta_init = NULL,
  mcp_max_iter = 5000L,
  mcp_tol = 1e-07,
  mcp_zero_tol = 1e-08,
  test_lambda = NULL,
  precision = c("nodewise", "diagonal", "user"),
  lambda_node = NULL,
  theta_rows = NULL,
  precision_fit = NULL,
  assume_independent = FALSE,
  alpha = 0.05,
  bootstrap_B = 999L,
  bootstrap_seed = NULL,
  bootstrap_multipliers = NULL,
  bootstrap_block_size = 256L,
  pilot_beta_init = NULL,
  pilot_operator_norm = NULL,
  pilot_max_iter = 30000L,
  pilot_tol = 5e-06,
  pilot_check_every = 50L,
  epsilon = 1e-10
)
```

## Arguments

- X:

  Numeric `n` by `p` design matrix.

- y:

  Numeric response vector.

- tested_set:

  One-based coordinates in the null hypothesis.

- t:

  Coordinatewise null target; defaults to zero.

- ridge_lambda:

  Optional fixed Ridge penalty.

- ridge_lambda_grid:

  Positive Ridge grid selected by GCV when `ridge_lambda` is `NULL`. GCV
  counts the unpenalized intercept in its effective degrees of freedom.

- mcp_lambda:

  Optional fixed MCP penalty.

- mcp_lambda_grid:

  Decreasing MCP grid selected by EBIC when `mcp_lambda` is `NULL`.

- mcp_gamma:

  MCP concavity parameter; the default is 3.

- mcp_ebic_gamma:

  EBIC combinatorial-penalty multiplier.

- mcp_nlambda:

  Number of automatically generated MCP penalties.

- mcp_lambda_min_ratio:

  Smallest automatic MCP penalty relative to the zero-solution penalty.

- mcp_beta_init:

  Optional original-scale MCP warm start. It never encodes or accepts an
  active set.

- mcp_max_iter:

  Maximum coordinate-descent sweeps per path point.

- mcp_tol:

  MCP convergence tolerance.

- mcp_zero_tol:

  Numerical support threshold.

- test_lambda:

  Square-root-LASSO penalty used only for the test pilot.

- precision:

  Precision-direction method for the test pilot.

- lambda_node:

  Nodewise-LASSO penalty when applicable.

- theta_rows:

  Optional user precision rows.

- precision_fit:

  Optional reusable precision fit for `tested_set`.

- assume_independent:

  Required for diagonal precision.

- alpha:

  Test level.

- bootstrap_B:

  Number of multiplier-bootstrap draws.

- bootstrap_seed:

  Optional bootstrap seed.

- bootstrap_multipliers:

  Optional reusable multiplier matrix.

- bootstrap_block_size:

  Bootstrap block size.

- pilot_beta_init:

  Optional original-scale square-root-LASSO warm start.

- pilot_operator_norm:

  Optional reusable design operator norm.

- pilot_max_iter:

  Maximum pilot iterations.

- pilot_tol:

  Pilot convergence tolerance.

- pilot_check_every:

  Pilot diagnostic interval.

- epsilon:

  Numerical floor in the Stein weight.

## Value

An `hd_shrinkage_fit` object whose `FM` endpoint is Ridge and whose `SM`
endpoint is the data-adaptive all-`X` MCP estimate. `PPS` is reported
separately from `PS`: it equals `FM` after rejection and `PS` otherwise.

## Details

With centered, standardized data, the full endpoint is
\$\$\hat\beta^{FM}=X^T(XX^T+n\lambda_R I_n)^{-1}y.\$\$ The sparse
endpoint minimizes \$\$\\y-Xb\\\_2^2/(2n)+\sum_j
p\_{\lambda_M,\gamma}^{MCP}(\|b_j\|).\$\$ The test pilot is deliberately
separate from both risk endpoints so that Ridge bias is not used to
calibrate the high-dimensional null test. Max-Stein risk dominance is
not asserted.

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(40 * 80), 40, 80)
beta <- c(1.2, -1, 0.8, rep(0, 77))
y <- drop(X %*% beta + rnorm(40))
fit <- fit_ridge_mcp_shrinkage(
  X, y, tested_set = 4:80,
  precision = "diagonal", assume_independent = TRUE,
  bootstrap_B = 99, mcp_nlambda = 20,
  pilot_max_iter = 5000, pilot_tol = 1e-4
)
coef(fit, method = "PS")
coef(fit, method = "PPS")
```
