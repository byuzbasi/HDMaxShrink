# Fit the CPSS-guided exact-null submodel endpoint

The selected core is refitted from scratch using an economy-SVD
Moore–Penrose operator. All coefficients outside the core are set
exactly to zero. The function never copies full-model coefficients and
never forms an OLS Gram inverse.

## Usage

``` r
fit_cpss_sm(X, y, core_set = NULL, selection = NULL, rank_tolerance = 1e-10)
```

## Arguments

- X:

  Numeric analysis-sample design matrix.

- y:

  Numeric analysis-sample response.

- core_set:

  One-based selected core indices, or a `cpss_core_selection` object.

- selection:

  Alternatively, an object returned by
  [`cpss_select_core()`](https://byuzbasi.github.io/HDMaxShrink/reference/cpss_select_core.md).
  Named features are aligned to `X` by identity.

- rank_tolerance:

  Relative economy-SVD rank tolerance.

## Value

A `cpss_sm_fit` object with original-scale coefficients and exact-null
diagnostics.

## References

Belloni, A. and Chernozhukov, V. (2013). Least squares after model
selection in high-dimensional sparse models. *Bernoulli*, 19, 521–547.
[doi:10.3150/11-BEJ410](https://doi.org/10.3150/11-BEJ410)
