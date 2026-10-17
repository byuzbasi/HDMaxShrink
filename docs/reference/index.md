# Package index

## Current CPSS–SM / Ridge workflow

- [`cpss_select_core()`](https://byuzbasi.github.io/HDMaxShrink/reference/cpss_select_core.md)
  : Select a stable high-dimensional submodel core with CPSS
- [`fit_cpss_sm()`](https://byuzbasi.github.io/HDMaxShrink/reference/fit_cpss_sm.md)
  : Fit the CPSS-guided exact-null submodel endpoint
- [`fit_cpss_ridge_shrinkage()`](https://byuzbasi.github.io/HDMaxShrink/reference/fit_cpss_ridge_shrinkage.md)
  : CPSS–SM shrinkage toward an all-predictor dual-Ridge endpoint

## Estimates and predictions

- [`coef(`*`<hd_shrinkage_fit>`*`)`](https://byuzbasi.github.io/HDMaxShrink/reference/coef.hd_shrinkage_fit.md)
  : Extract coefficients from a high-dimensional shrinkage fit
- [`predict(`*`<hd_shrinkage_fit>`*`)`](https://byuzbasi.github.io/HDMaxShrink/reference/predict.hd_shrinkage_fit.md)
  : Predict from a high-dimensional shrinkage fit

## Retained frameworks and inference

- [`fit_hd_shrinkage()`](https://byuzbasi.github.io/HDMaxShrink/reference/fit_hd_shrinkage.md)
  : Restriction-adaptive shrinkage in a high-dimensional linear model
- [`fit_partial_sqrt_lasso_shrinkage()`](https://byuzbasi.github.io/HDMaxShrink/reference/fit_partial_sqrt_lasso_shrinkage.md)
  : High-dimensional square-root-LASSO shrinkage toward an exact
  coordinate null
- [`fit_ridge_mcp_shrinkage()`](https://byuzbasi.github.io/HDMaxShrink/reference/fit_ridge_mcp_shrinkage.md)
  : Tmax-guided shrinkage between all-X Ridge and all-X MCP estimators
- [`sqrt_lasso_hd()`](https://byuzbasi.github.io/HDMaxShrink/reference/sqrt_lasso_hd.md)
  : Fit a square-root lasso model with an RcppArmadillo solver
- [`estimate_precision_rows()`](https://byuzbasi.github.io/HDMaxShrink/reference/estimate_precision_rows.md)
  : Estimate target rows of a high-dimensional precision matrix
- [`max_score_test_hd()`](https://byuzbasi.github.io/HDMaxShrink/reference/max_score_test_hd.md)
  : Conditional-Gaussian exact-null maximum partial-t test
- [`max_test_hd()`](https://byuzbasi.github.io/HDMaxShrink/reference/max_test_hd.md)
  : Multiplier-bootstrap maximum test from debiased influence scores

## Package and low-level numerical utilities

- [`HDMaxShrink`](https://byuzbasi.github.io/HDMaxShrink/reference/HDMaxShrink-package.md)
  [`HDMaxShrink-package`](https://byuzbasi.github.io/HDMaxShrink/reference/HDMaxShrink-package.md)
  : HDMaxShrink: test-guided high-dimensional shrinkage
