# Data, reproducibility and validation

## Input dictionary

| Input | Meaning | Required separation |
|:---|:---|:---|
| `X_selection`, `y_selection` | Numeric predictor matrix and response used to select a core | Independent of analysis observations |
| `X_analysis`, `y_analysis` | Numeric design and response used to estimate and test | Matching predictor names and order |
| Mandatory core | Predictor identities fixed in advance | Not selected using analysis responses |
| Test design and response | Untouched evaluation observations | Never used to choose tuning parameters |
| `bootstrap_seed` | Reproducible conditional-calibration seed | Record it with the draw count and version |

Rows are observations; columns are predictors. Existing fit routines
perform their documented centering and RMS scaling internally and return
coefficients on the original scale. Do not add undocumented cross-sample
preprocessing.

## Reproduce an example, not a production claim

The [quick-start
script](https://github.com/byuzbasi/HDMaxShrink/blob/main/inst/examples/quickstart.R)
and saved figure belong to one simulated dataset. The package also
retains historical benchmark sources under
[inst/benchmarks](https://github.com/byuzbasi/HDMaxShrink/tree/main/inst/benchmarks).
Their presence is not a statement that every study was rerun for this
site.

For an analysis record, save the call, seeds, package version,
data-source identifiers, selection/analysis/test split and
convergence/calibration output. Record
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) alongside
results. Do not overwrite historical outputs.

## Citation and contribution

Use `citation("HDMaxShrink")` for software and methods records, and cite
the preprint as a preprint. Report bugs with a small synthetic example
and [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html); do not
attach private study data. The website build renders documentation and
existing figures. It does not run production calibration, Monte Carlo
studies or performance benchmarks.
