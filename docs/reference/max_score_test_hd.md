# Conditional-Gaussian exact-null maximum partial-t test

This convenience function runs exactly the same full-model pilot and
maximum test used by
[`fit_partial_sqrt_lasso_shrinkage()`](https://byuzbasi.github.io/HDMaxShrink/reference/fit_partial_sqrt_lasso_shrinkage.md)
and returns the inference component. The default is the
conditional-Gaussian exact-null maximum partial-t test; it is not a
Wald, chi-square, or joint partial-F test.

## Usage

``` r
max_score_test_hd(
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

A maximum-test list with the fitted pilot attached as `full_model_beta`
and `full_solver`.
