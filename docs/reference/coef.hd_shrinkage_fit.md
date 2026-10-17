# Extract coefficients from a high-dimensional shrinkage fit

Extract coefficients from a high-dimensional shrinkage fit

## Usage

``` r
# S3 method for class 'hd_shrinkage_fit'
coef(
  object,
  method = c("PS", "PPS", "FM", "SM", "PT", "S", "RIDGE", "MCP", "FSL", "SMSL", "FDSL",
    "SMDSL", "full", "submodel", "restricted"),
  ...
)
```

## Arguments

- object:

  A fitted `hd_shrinkage_fit` object.

- method:

  Selected endpoint (`"FM"` or `"SM"`), adaptive estimator (`"PT"`,
  `"S"`, `"PS"`, or `"PPS"`), raw square-root lasso endpoint (`"FSL"` or
  `"SMSL"`), explicit debiased endpoint (`"FDSL"` or `"SMDSL"`) when
  available, or the all-`X` endpoints `"RIDGE"` and `"MCP"` from
  [`fit_ridge_mcp_shrinkage()`](https://byuzbasi.github.io/HDMaxShrink/reference/fit_ridge_mcp_shrinkage.md).
  The legacy names `"full"` and `"restricted"` remain available.

- ...:

  Unused.

## Value

A named coefficient vector including the intercept.
