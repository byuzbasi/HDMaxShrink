# HDMaxShrink 0.6.2

- Replaced the optional CPSS--MCP base-path dependency on `ncvreg` by
  singleton-group `grpreg` MCP. With one predictor per group, the group MCP
  penalty is the ordinary coordinatewise MCP penalty.
- Corrected `lambda_min_ratio` provenance: both CPSS--LASSO and CPSS--MCP now
  record a lower path endpoint relative to the engine-specific maximum
  penalty. The earlier absolute-`lambda.min` label was documentation metadata,
  not the behavior of the MCP engine.
- Retained the `selected_path_locally_convex` diagnostic as `NA` for backward
  compatibility because `grpreg` does not expose the former local-convexity
  index. No local-convexity claim is made in version 0.6.2.

# HDMaxShrink 0.6.1

- Added portable package-level LAPACK, BLAS, and Fortran linkage through
  `src/Makevars` and `src/Makevars.win`. This fixes unresolved Armadillo
  eigendecomposition and SVD symbols on Linux without changing estimators,
  APIs, defaults, seeds, or statistical calculations.
- Extended source-provenance checks to include the package build recipes.

# HDMaxShrink 0.6.0

- Added backward-compatible `mandatory_core` and `candidate_set` arguments to
  `cpss_select_core()`. In every complementary half-sample, required
  predictors are centered and projected out by a checked economy SVD before
  LASSO or MCP sees the optional candidate family.
- Separated mandatory variables, CPSS-selected extensions, ineligible
  optional directions, and predictors outside the selection universe in the
  returned audit trail. Empty extensions now leave a nonempty mandatory core
  intact without a top-k fallback.
- Scoped the reported Meinshausen--Buehlmann expression to the eligible
  optional family; mandatory variables have no artificial CPSS frequency.
- Updated the LASSO base fit for the `glmnet` 5.0 control interface while
  retaining compatibility with earlier interfaces.

# HDMaxShrink 0.5.0

- Added `cpss_select_core()` for strict complementary-pairs stability
  selection with a size-budgeted LASSO base selector and an optional MCP
  sensitivity selector. The implementation records feature identities,
  stability frequencies, base-fit convergence, MCP local-convexity
  diagnostics, and the applicability of the above-half stability-selection
  PFER expression. It never silently replaces an empty stable core by top-k
  variables.
- Added `fit_cpss_sm()`, a CPSS-guided exact-null submodel estimator computed
  by an economy-SVD Moore--Penrose refit on an analysis sample independent of
  selection.
- Added `fit_cpss_ridge_shrinkage()`, coupling CPSS--SM with an all-predictor
  dual-Ridge full-model endpoint and the existing conditional-Gaussian
  max-partial-t/inverse-moment PT, S, PS, and PPS family. Ridge is not used to
  calibrate the null test, and no OLS, p-by-p, or q-by-q inverse is formed.

# HDMaxShrink 0.4.0

- Added `full_endpoint = "all_x_square_root_lasso"` as the primary endpoint.
  It centers and RMS-scales all columns of `X` and `y`, fixes the standardized
  intercept at zero, and estimates all `p` coefficients jointly in one
  square-root-LASSO optimization with positive penalties on every coordinate.
- Kept `SM` as a fresh exact-null SVD refit using only the prespecified
  low-dimensional `X_A`; no FM coefficient is copied into SM.
- Back-transforms FM, SM, PT, S, PS, and PPS coefficients and predictions to
  the original data scale while recording the zero standardized intercept.
- Preserved the conditional-Gaussian maximum partial-t test on the
  `X_A`-residualized tested block. The all-X fit does not change the null-test
  statistic or its Monte Carlo calibration.
- Retained the thresholded partial-debiased and profiled partial endpoints as
  explicit backward-compatible research modes; they are no longer the
  primary analysis endpoint.

# HDMaxShrink 0.3.0

- Added `fit_partial_sqrt_lasso_shrinkage()` as the primary non-Ridge
  framework: a partial square-root LASSO supplies the pilot, `FM` is its
  thresholded partial-debiased and core-reprofiled endpoint, and `SM` is the
  exact-null SVD refit under `beta[B] = 0`.
- Added an explicit marginal exceedance budget `threshold_eta`; the
  `q = 990` study prespecifies eta 5 and reports eta 1, 3, and 5 sensitivity.
- Profiled the prespecified low-dimensional core by an economy SVD. No OLS
  Gram inverse, `q` by `q` inverse, or `p` by `p` inverse is formed.
- Added a Gaussian-iid diagonal one-step `Tmax` and multiplier-bootstrap null
  calibration, exposed separately through `max_score_test_hd()`.
- Synchronized the primary partial fit with the exact-null Gaussian theory:
  its default is now the conditional-Gaussian maximum partial-t test and the
  reciprocal inverse-moment max-Stein calibration. The raw score maximum,
  partial-t maximum, null second moment, and inverse-moment calibration are
  stored in separate fields. The former residual-multiplier plus second-moment
  construction remains available through explicit arguments.
- Added reusable null geometry, full-model warm starts, and tests showing that
  the null endpoint is a genuine refit rather than a projection or oracle
  construction.
- Retained Ridge--MCP and the earlier projection framework without changing
  their APIs.

# HDMaxShrink 0.2.0

- Added a dense--sparse endpoint framework with all-`X` dual Ridge as `FM` and
  all-`X` MCP as `SM`; the MCP solver receives no active-set information.
- Added RcppArmadillo Ridge-path and MCP coordinate-descent kernels, with GCV
  and EBIC tuning respectively.
- Kept the large-block test on a separate debiased square-root-LASSO score so
  that Ridge bias does not calibrate `Tmax`.
- Retained the 0.1.2 projection estimator as a backward-compatible research
  mode.

# HDMaxShrink 0.1.2

- Replaced the reduced-design submodel fit with an exact restriction
  projection of the single all-`X` full square-root lasso estimator.
- Projected the all-`X` debiased full endpoint instead of separately
  debiasing a reduced submodel.
- Kept the multiplier-bootstrap `Tmax` construction unchanged and added
  projection-invariance tests.

# HDMaxShrink 0.1.1

- Added symmetric one-step debiasing of both full-model and submodel fits.
- Added explicit FM, SM, FSL, SMSL, FDSL, and SMDSL coefficient selectors.
- Preserved the regularized endpoints as the backward-compatible default.

# HDMaxShrink 0.1.0

- Added RcppArmadillo primal-dual square-root lasso kernels.
- Added a memory-efficient zero-block restricted solver.
- Added small-`q` debiased Wald inference.
- Added large-block multiplier-bootstrap maximum inference for `q > n`.
- Added preliminary-test, Stein-type, and positive-part max-Stein estimators.
