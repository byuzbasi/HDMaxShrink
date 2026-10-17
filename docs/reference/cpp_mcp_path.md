# MCP coordinate-descent path on all design columns

The objective is `||y-X beta||^2/(2n)` plus the coordinatewise MCP
penalty. Successive decreasing penalties use warm starts.

## Usage

``` r
cpp_mcp_path(
  X,
  y,
  lambda,
  beta_init,
  gamma = 3,
  max_iter = 5000L,
  tol = 1e-07,
  zero_tol = 1e-08
)
```
