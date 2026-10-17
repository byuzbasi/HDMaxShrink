# Apply cached core-model geometry to a centered response

The returned influence scores are the heteroskedasticity-robust,
residualized score contributions used by the multiplier-bootstrap
maximum test of the coordinate null.

## Usage

``` r
cpp_partial_null_apply(
  y,
  coefficient_operator,
  orthonormal_basis,
  residualized_tested
)
```
