# Restriction-adaptive shrinkage in a high-dimensional linear model

`fit_hd_shrinkage()` fits one full square-root lasso with all columns of
`X`. The submodel endpoint is the restriction projection of that full
estimator; it is never re-fitted on a reduced design. The endpoints can
be the regularized full/projection pair or, for a zero-coordinate
restriction, a full one-step debiased estimator and its projection. A
debiased Wald statistic is available for fixed, small restriction
dimension. A multiplier-bootstrap maximum test is used for a large
coordinate block and does not form a `q` by `q` covariance matrix.
Neither branch uses ordinary, ridge, or generalized ridge regression.

## Usage

``` r
fit_hd_shrinkage(
  X,
  y,
  restricted_set,
  C = NULL,
  t = NULL,
  endpoint = c("regularized", "debiased_models"),
  lambda = NULL,
  lambda_node = NULL,
  test = c("auto", "wald", "max"),
  precision = c("nodewise", "diagonal", "user"),
  theta_rows = NULL,
  precision_fit = NULL,
  full_model_precision_fit = NULL,
  submodel_precision_fit = NULL,
  assume_independent = FALSE,
  alpha = 0.05,
  bootstrap_B = 999L,
  bootstrap_seed = NULL,
  bootstrap_multipliers = NULL,
  bootstrap_block_size = 256L,
  fixed_q_max = 10L,
  penalty_factor = NULL,
  max_iter = 30000L,
  tol = 5e-06,
  check_every = 50L,
  verbose = FALSE,
  full_beta_init = NULL,
  restricted_beta_init = NULL,
  full_operator_norm = NULL,
  restricted_operator_norm = NULL,
  epsilon = 1e-10
)
```

## Arguments

- X:

  Numeric `n` by `p` design matrix.

- y:

  Numeric response vector of length `n`.

- restricted_set:

  One-based indices `M` involved in the restriction.

- C:

  Optional `q` by `length(M)` contrast matrix. `NULL` means the
  coordinate restriction `beta[M] = t`, whose projection is implemented
  without constructing a `q` by `q` matrix.

- t:

  Null target. It defaults to zero.

- endpoint:

  Shrinkage endpoints. `"regularized"` uses the full square-root lasso
  and its restriction projection. `"debiased_models"` computes one
  full-model debiased estimator and projects it to the coordinate zero
  subspace. The latter currently requires a coordinate zero-block
  restriction.

- lambda:

  Square-root lasso penalty.

- lambda_node:

  Nodewise lasso penalty.

- test:

  One of `"auto"`, `"wald"`, or `"max"`. Automatic selection uses Wald
  only when `q <= fixed_q_max`.

- precision:

  One of `"nodewise"`, `"diagonal"`, or `"user"`.

- theta_rows:

  User precision rows for `precision = "user"`.

- precision_fit:

  Optional reusable result from
  [`estimate_precision_rows()`](https://byuzbasi.github.io/HDMaxShrink/reference/estimate_precision_rows.md).

- full_model_precision_fit:

  Optional reusable all-row precision result for the standardized full
  design when `endpoint = "debiased_models"`.

- submodel_precision_fit:

  Deprecated compatibility argument. It is ignored because no separate
  submodel debiasing fit is performed.

- assume_independent:

  Required for diagonal precision. This is suitable for a design whose
  diagonal population precision is assumed in advance, such as the
  Gaussian iid simulation design.

- alpha:

  Test level.

- bootstrap_B:

  Number of multiplier draws for the maximum test.

- bootstrap_seed:

  Optional seed used without changing the caller's RNG state.

- bootstrap_multipliers:

  Optional reusable `n` by `bootstrap_B` Gaussian multiplier matrix.

- bootstrap_block_size:

  Number of multiplier draws processed in one C++ block.

- fixed_q_max:

  Largest restriction dimension selected for Wald by `test = "auto"`.

- penalty_factor:

  Nonnegative square-root lasso penalty factors.

- max_iter:

  Maximum primal-dual iterations.

- tol:

  Solver tolerance.

- check_every:

  Iteration interval for convergence diagnostics.

- verbose:

  Print periodic solver diagnostics.

- full_beta_init:

  Optional original-scale warm start for the full fit.

- restricted_beta_init:

  Deprecated compatibility argument. It is ignored because the
  restricted endpoint is projected from the full fit.

- full_operator_norm:

  Optional reusable full-design operator norm.

- restricted_operator_norm:

  Deprecated compatibility argument. It is ignored because no
  restricted-design solver is run.

- epsilon:

  Positive numerical floor used only in a shrinkage denominator.

## Value

An object of class `hd_shrinkage_fit` containing the regularized and
selected full/submodel endpoints, the available adaptive estimators,
full solver and projection diagnostics, precision directions, and
inference output.

## Details

Let `A beta = t` denote the restriction and let `V` denote the
projection geometry. The conceptual restricted endpoint is
\$\$\hat\beta^{SM}=\hat\beta^F- V A^T(A V
A^T)^\dagger(A\hat\beta^F-t).\$\$ The implementation uses Euclidean
geometry after internal column standardization. For `beta[M] = t`, this
simply replaces the selected coordinates of the all-`X` full estimate by
`t`; all other coordinates are inherited exactly from the full fit.

For a large coordinate restriction, the max-calibrated estimators are
\$\$\hat\beta^{S\_\infty}=\hat\beta^{SM}+
(1-a^\*/(T\_{\max}^2\vee\epsilon)) (\hat\beta^{FM}-\hat\beta^{SM})\$\$
and its positive-part version. The separately reported protected version
equals the full endpoint when the max test rejects and otherwise equals
the positive-part estimator. Here `a*` is the conditional mean of the
squared multiplier-bootstrap maximum. These estimators are experimental;
the package does not assert classical James–Stein dominance for them.

The restricted set must be specified independently of the response used
for inference. If it is selected from data, sample splitting or external
information is required for the advertised test interpretation.

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(40 * 70), 40, 70)
beta <- c(rep(1, 4), rep(0, 66))
y <- drop(X %*% beta + rnorm(40))
fit <- fit_hd_shrinkage(
  X, y,
  restricted_set = 5:70,
  test = "max",
  precision = "diagonal",
  assume_independent = TRUE,
  bootstrap_B = 99,
  max_iter = 5000,
  tol = 1e-4
)
coef(fit, method = "PS")
coef(fit, method = "PPS")
```
