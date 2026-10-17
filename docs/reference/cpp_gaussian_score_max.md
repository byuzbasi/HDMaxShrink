# Conditional Gaussian score maxima for an exact coordinate null

The input residual must be the centered exact-null residual, the tested
design must already be residualized off the centered core, and the
columns of `core_basis` must be an orthonormal basis for that core. Each
Gaussian draw is centered and projected off the core before applying the
same residual-scale studentization as the observed statistic.

## Usage

``` r
cpp_gaussian_score_max(
  null_residual,
  residualized_tested,
  core_basis,
  gaussian_draws,
  block_size = 256L
)
```
