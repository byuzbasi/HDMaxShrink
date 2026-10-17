# Select a stable high-dimensional submodel core with CPSS

`cpss_select_core()` implements complementary-pairs stability selection
(CPSS) with a size-budgeted LASSO or MCP base selector. Selection must
be performed on data independent of any sample later used for the
conditional-Gaussian exact-null test. CPSS itself is due to Shah and
Samworth (2013); this function supplies the core-selection stage used by
the package's CPSS–SM endpoint.

## Usage

``` r
cpss_select_core(
  X,
  y,
  selector = c("lasso", "mcp"),
  complementary_pairs = 50L,
  base_selection_size = 20L,
  stability_threshold = 0.6,
  seed = NULL,
  path_points = 100L,
  lambda_min_ratio = 0.01,
  coefficient_tolerance = 1e-08,
  lasso_maxit = 100000L,
  mcp_gamma = 3,
  mcp_max_iter = 100000L,
  no_top_k_fallback = TRUE,
  strata = NULL,
  mandatory_core = NULL,
  candidate_set = NULL,
  rank_tolerance = 1e-10
)
```

## Arguments

- X:

  Numeric selection-sample design matrix.

- y:

  Numeric selection-sample response.

- selector:

  Base selector, either `"lasso"` or `"mcp"`.

- complementary_pairs:

  Number of complementary half-sample pairs.

- base_selection_size:

  Maximum support size of each base fit.

- stability_threshold:

  Inclusion-frequency threshold. The usual Meinshausen–Buehlmann PFER
  expression is only reported above one half.

- seed:

  Optional reproducible resampling seed.

- path_points:

  Number of penalty-path points per base fit.

- lambda_min_ratio:

  Path lower-bound control. It is passed as `lambda.min.ratio` to
  `glmnet` for LASSO and as the `lambda.min` fraction of `lambda.max` to
  `grpreg` for the singleton-group MCP sensitivity path.

- coefficient_tolerance:

  Numerical support threshold.

- lasso_maxit:

  Maximum `glmnet` iterations.

- mcp_gamma:

  MCP concavity parameter.

- mcp_max_iter:

  Maximum total `grpreg` iterations over the MCP path.

- no_top_k_fallback:

  Must be `TRUE`. An empty thresholded core is returned without a
  data-dependent top-k replacement. If a mandatory core is supplied, an
  empty CPSS extension leaves the mandatory core unchanged.

- strata:

  Optional resampling strata with one entry per row.

- mandatory_core:

  Optional one-based indices or feature names that are always retained.
  Within every half-sample they are centered and projected out by a
  checked economy SVD before the base selector is fitted. Mandatory
  variables are outside the CPSS frequency and PFER family.

- candidate_set:

  Optional one-based indices or feature names eligible for CPSS
  selection. The default is every predictor outside `mandatory_core`.

- rank_tolerance:

  Relative SVD tolerance for mandatory-core residualization.

## Value

A `cpss_core_selection` object containing the selected indices,
selection frequencies, feature identities, and base-fit diagnostics.

## References

Meinshausen, N. and Buehlmann, P. (2010). Stability selection. *Journal
of the Royal Statistical Society: Series B*, 72, 417–473.
[doi:10.1111/j.1467-9868.2010.00740.x](https://doi.org/10.1111/j.1467-9868.2010.00740.x)

Shah, R. D. and Samworth, R. J. (2013). Variable selection with error
control: another look at stability selection. *Journal of the Royal
Statistical Society: Series B*, 75, 55–80.
[doi:10.1111/j.1467-9868.2011.01034.x](https://doi.org/10.1111/j.1467-9868.2011.01034.x)
