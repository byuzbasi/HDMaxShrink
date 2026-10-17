# High-dimensional square-root-LASSO shrinkage toward an exact coordinate null

This is a non-Ridge framework for a sparse linear model with `p > n`.
The columns are split in advance into a low-dimensional core `A` and a
high-dimensional tested block `B`. The recommended all-X endpoint
centers and scales both `X` and `y`, fixes the standardized intercept at
zero, and estimates every coefficient jointly in one square-root-LASSO
problem. Legacy partial endpoints remain available for reproducibility.
The submodel endpoint is the exact-null fit under `beta[B] = 0`; it is
not an oracle fit and is not obtained by truncating the full estimate.

## Usage

``` r
fit_partial_sqrt_lasso_shrinkage(
  X,
  y,
  core_set,
  tested_set = NULL,
  lambda = NULL,
  penalty_loadings = NULL,
  all_x_penalty_loadings = NULL,
  standardize_y = NULL,
  threshold_eta = NULL,
  assume_independent = FALSE,
  alpha = 0.05,
  bootstrap_B = 999L,
  bootstrap_seed = NULL,
  bootstrap_multipliers = NULL,
  bootstrap_block_size = 256L,
  max_iter = 30000L,
  tol = 5e-06,
  check_every = 50L,
  verbose = FALSE,
  full_beta_init = NULL,
  full_operator_norm = NULL,
  null_geometry = NULL,
  rank_tolerance = 1e-10,
  epsilon = 1e-10,
  test_calibration = c("conditional_gaussian", "residual_multiplier"),
  shrinkage_calibration = c("inverse_moment", "second_moment"),
  gaussian_draws = NULL,
  full_endpoint = c("thresholded_partial_debiased", "profiled_square_root_lasso",
    "all_x_square_root_lasso")
)
```

## Arguments

- X:

  Numeric `n` by `p` design matrix.

- y:

  Numeric response vector.

- core_set:

  Prespecified one-based indices of the low-dimensional core. It must
  contain fewer than `n` columns and must not be chosen from `y` without
  sample splitting.

- tested_set:

  Prespecified high-dimensional null block. By default it is the
  complement of `core_set`; together the two sets must partition `X`.

- lambda:

  Positive square-root-LASSO penalty. The default is
  `1.1 * sqrt(2 * log(2 * p) / n)` for the all-X endpoint and
  `1.1 * sqrt(2 * log(2 * q) / n)` for a legacy partial endpoint.

- penalty_loadings:

  Optional positive loadings for the tested block. The default is the
  empirical norm of each tested column after projection off the core.

- all_x_penalty_loadings:

  Optional length-`p` positive loadings for
  `full_endpoint = "all_x_square_root_lasso"`. The default penalizes all
  columns equally. No core coefficient is exempted from the FM penalty.

- standardize_y:

  Whether to center and scale `y` before fitting. It is fixed at `TRUE`
  for the all-X endpoint and must remain `FALSE` for the legacy
  endpoints. `X` is always centered and scaled. The standardized
  intercept is exactly zero; returned coefficients and prediction
  intercepts are back-transformed to the original data scale.

- threshold_eta:

  Positive marginal exceedance budget for the thresholded debiased full
  endpoint. The threshold is `qnorm(1 - threshold_eta / (2 * q))`. The
  default is `min(5, q / 2)`, which gives the prespecified finite-sample
  budget `5` in the advertised ultra-high-dimensional designs. This is a
  PFER-style tuning interpretation, not an FDR guarantee. The vanishing
  sequence `1 / log(q)` is retained for theoretical and sensitivity
  analyses. This argument must remain `NULL` for either unthresholded
  square-root-LASSO endpoint because no coefficient threshold is applied
  there.

- assume_independent:

  Must be `TRUE` for the default thresholded partial-debiased FM
  endpoint. It records the Gaussian-iid (or otherwise prespecified
  diagonal population precision) condition used by that endpoint. It is
  not required by either unthresholded square-root-LASSO endpoint or by
  the conditional-Gaussian exact-null test.

- alpha:

  Test level.

- bootstrap_B:

  Number of Gaussian multiplier draws.

- bootstrap_seed:

  Optional seed that leaves the caller RNG unchanged.

- bootstrap_multipliers:

  Optional reusable `n` by `bootstrap_B` matrix. Under the default
  conditional-Gaussian calibration it is accepted as a
  backward-compatible alias for `gaussian_draws`; under
  `test_calibration = "residual_multiplier"` it contains the multiplier
  draws used by the legacy test.

- bootstrap_block_size:

  Number of draws processed per C++ block.

- max_iter:

  Maximum primal-dual iterations.

- tol:

  Solver convergence tolerance.

- check_every:

  Iteration interval for diagnostics.

- verbose:

  Print periodic solver diagnostics.

- full_beta_init:

  Optional original-scale length-`p` warm start.

- full_operator_norm:

  Optional cached operator norm of the full standardized design for the
  all-X endpoint or of the residualized tested design for a legacy
  partial endpoint.

- null_geometry:

  Optional cached geometry returned in a previous fit on exactly the
  same design and partition.

- rank_tolerance:

  Relative numerical-rank tolerance for the core SVD.

- epsilon:

  Positive numerical floor in the Stein denominator.

- test_calibration:

  Test calibration. The default `"conditional_gaussian"` is the
  finite-sample exact-null Gaussian max-partial-t calibration.
  `"residual_multiplier"` retains the legacy partial-debiased multiplier
  maximum.

- shrinkage_calibration:

  Max-Stein calibration. The default `"inverse_moment"` uses the
  reciprocal bootstrap inverse second moment. `"second_moment"` retains
  the legacy bootstrap second-moment rule.

- gaussian_draws:

  Optional reusable `n` by `bootstrap_B` matrix for the
  conditional-Gaussian exact-null calibration.

- full_endpoint:

  Definition of the full-model endpoint. The default
  `"thresholded_partial_debiased"` preserves the Gaussian-iid framework.
  `"profiled_square_root_lasso"` uses the reprofiled partial
  square-root-LASSO pilot directly and is available for correlated fixed
  designs. `"all_x_square_root_lasso"` estimates all `p` coefficients in
  one centered-and-scaled square-root-LASSO problem and is the
  recommended endpoint when FM must use the complete design without an
  A/B-specific estimation step. Both unthresholded endpoints require
  `test_calibration = "conditional_gaussian"`.

## Value

An object of class `hd_shrinkage_fit`. `FM` is the selected full
endpoint, `SM` is the exact-null endpoint, and `PT`, `S`, `PS`, and
`PPS` are test-guided estimators between those endpoints. `PPS` equals
`FM` after rejection and otherwise equals `PS`; it is reported as a
separate protected estimator and does not replace `PS`.
`full_regularized` equals `FM` for either unthresholded
square-root-LASSO endpoint.

## Details

For the recommended all-X endpoint, write `X*` and `y*` for columnwise
and responsewise centered-and-RMS-scaled training data. It solves once
\$\$\widehat\beta^{FM,\*}=\arg\min\_{b\in\mathbb R^p}
\\y^\*-X^\*b\\\_2/\sqrt n+\lambda\sum\_{j=1}^p\omega_j\|b_j\|.\$\$ Its
standardized intercept is exactly zero. Coefficients and the prediction
intercept are then mapped back to the original scale. The A/B partition
is used only to construct the exact-null SM and the restriction test, so
changing that partition does not change the all-X FM for fixed `X`, `y`,
penalty, and solver settings.

For a legacy partial endpoint, let `M_A` project orthogonally off the
core, `r_A = M_A y`, and `V = M_A X_B`. The pilot solves \$\$\hat
b_B^{FM}=\arg\min_b \\r_A-Vb\\\_2/\sqrt
n+\lambda\sum_j\omega_j\|b_j\|,\$\$ and produces one-step coordinates
`tilde beta_j` and standard errors. Let `d = n - rank(X_A) - 1` be the
residual degrees of freedom after centering and projecting off the core.
With \$\$Z_j=\sqrt d\\\widetilde\beta_j/\widehat{se}\_j,\qquad
\tau\_\eta=\Phi^{-1}(1-\eta/(2q)),\$\$ the full endpoint is
\$\$\widehat\beta\_{B,j}^{FM}=\widetilde\beta_j
1\\\|Z_j\|\>\tau\_\eta\\,\qquad
\widehat\beta_A^{FM}=X_A^+(y-X_B\widehat\beta_B^{FM}).\$\$ Under
`full_endpoint = "profiled_square_root_lasso"`, the tested block is
instead the optimizer `hat b_B` in the first display and
\$\$\widehat\beta_A^{FM}=X_A^+(y-X_B\hat b_B).\$\$ This alternative does
not debias or threshold the tested coefficients and therefore does not
invoke a diagonal precision approximation. The exact-null endpoint is
\$\$\hat b^{SM}=(X_A^+y,0_q).\$\$ The Moore–Penrose operator is computed
from a checked economy SVD; no OLS Gram inverse and no `p` by `p`
inverse is used.

Conditional on the fixed design, the default test starts from the
exact-null residual, applies the same core projection to every Gaussian
draw, and transforms the maximum normalized score monotonically to the
maximum partial-t scale. Its Monte Carlo p-value is therefore unchanged
by the transformation. No direct chi-square limit and no `q` by `q`
covariance inverse is asserted. By default the max-Stein factor is
`1 - kappa_star / T_max^2`, where
`kappa_star = 1 / mean(1 / (T_max_star)^2)`; the positive-part factor
truncates this at zero. With `R` denoting the rejection indicator, the
protected weight is `R + (1 - R) * positive_weight`. Classical
James–Stein dominance is not claimed for either estimator. Finite-sample
conditional exactness requires homoskedastic Gaussian errors and a
prespecified core, but it does not require independent design columns.

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(50 * 100), 50, 100)
beta <- c(1.2, -1, 0.8, rep(0, 97))
y <- drop(X %*% beta + rnorm(50))
fit <- fit_partial_sqrt_lasso_shrinkage(
  X, y, core_set = 1:3, tested_set = 4:100,
  assume_independent = TRUE, bootstrap_B = 99,
  max_iter = 5000, tol = 1e-4
)
coef(fit, method = "PS")
coef(fit, method = "PPS")
fit_all_x <- fit_partial_sqrt_lasso_shrinkage(
  X, y, core_set = 1:3, tested_set = 4:100,
  bootstrap_B = 99, max_iter = 5000, tol = 1e-4,
  full_endpoint = "all_x_square_root_lasso"
)
```
