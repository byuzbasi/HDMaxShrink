# Dual high-dimensional ridge path

Computes ridge estimates from an `n` by `n` eigendecomposition rather
than a `p` by `p` inverse. The objective is
`||y-X beta||^2/(2n) + lambda ||beta||^2/2`.

## Usage

``` r
cpp_ridge_dual_path(X, y, lambda)
```
