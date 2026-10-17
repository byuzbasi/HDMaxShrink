# SVD geometry for a prespecified low-dimensional core model

The routine never forms `(X_core' X_core)^{-1}`. It returns the
Moore–Penrose coefficient operator (after verifying full column rank),
an orthonormal basis for the core column space, and every tested column
residualized against that space.

## Usage

``` r
cpp_partial_null_geometry(X_core, X_tested, rank_tolerance = 1e-10)
```
