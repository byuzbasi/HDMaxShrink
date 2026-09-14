# HDMaxShrink

## Install and try

Source repository: <https://github.com/byuzbasi/HDMaxShrink>.

```r
# Install remotes first if it is not already available.
install.packages("remotes")
remotes::install_github("byuzbasi/HDMaxShrink", upgrade = "never")
```

R >= 4.1.0 and a working C++17 toolchain are required. The package links to
R's LAPACK/BLAS and uses RcppArmadillo. The optional CPSS-MCP selector also
requires `grpreg`, listed under Suggests; install it separately if needed.

Run the small, self-contained example after installation:

```r
source(system.file("examples", "quickstart.R", package = "HDMaxShrink"))
```

The example uses independent selection and analysis samples with 120
predictors and 60 observations per sample. It is an executable demonstration,
not a reproduction of the manuscript's production study. It returns FM, SM,
PT, S and PS coefficient estimates.

Report software problems at <https://github.com/byuzbasi/HDMaxShrink/issues>.
Use `citation("HDMaxShrink")` for the software reference. The package is
licensed under GPL version 3 or later, as declared in `DESCRIPTION`.

This public source snapshot retains the frozen 0.6.2 numerical implementation
and tests. GitHub metadata, installation guidance and a small example were
added for distribution; the original study's versioned source archives remain
unchanged. This repository contains the R package, not raw data, manuscript
submission files or the full study's execution checkpoints.

HDMaxShrink is a research package for pretest--Stein estimation in sparse
linear models with \(p\gg n\). Version 0.6.2 preserves the version 0.6.1
statistical API and uses singleton-group `grpreg` MCP for the optional
CPSS--MCP sensitivity path. Version 0.6.1 added portable linkage for the
compiled numerical kernels, and version 0.6.0 added a reusable mandatory-core
conditional CPSS workflow. CPSS and post-selection refitting are established
methods; the proposed contribution is their honest integration with an
exact-null SVD submodel, an all-\(X\) dual-Ridge full endpoint, and a
conditionally calibrated max-partial-\(t\) PT/Stein family.

## CPSS-guided submodel and Ridge full model

Use an independent selection sample to learn a stable optional extension
conditional on required predictors \(A_0\):

    selection <- cpss_select_core(
      X_selection, y_selection,
      selector = "lasso",       # "mcp" is a sensitivity selector
      complementary_pairs = 50,
      base_selection_size = 20,
      stability_threshold = 0.60,
      seed = 20260825,
      mandatory_core = c("lineage_2", "lineage_3", "MSIScore"),
      candidate_set = grep("^expr::", colnames(X_selection))
    )

Then use a separate analysis sample:

    fit <- fit_cpss_ridge_shrinkage(
      X_analysis, y_analysis,
      selection = selection,
      ridge_lambda = 0.25,
      bootstrap_B = 999,
      bootstrap_seed = 20260826
    )

    coef(fit, method = "FM")
    coef(fit, method = "SM")
    coef(fit, method = "PS")

Within every complementary half, the response and optional candidate columns
are centered and residualized against the mandatory design before the base
selector is fitted. Mandatory variables are outside the CPSS frequency and
PFER family. The final core is
\(\widehat A=A_0\cup\widehat E\); if the selected extension is empty it equals
\(A_0\), with no top-k fallback.

FM estimates all columns of \(X\) jointly by direct dual Ridge. SM is
recomputed only on \(X_{\widehat A}\) with an economy-SVD Moore--Penrose
operator and sets every complementary coefficient exactly to zero. The null
test uses the exact-null SM residual; Ridge bias never calibrates the test.
Conditional finite-sample exactness requires that selection be independent of
the analysis response and that the selected null is true, in addition to iid
homoskedastic Gaussian errors. CPSS--MCP is reported as sensitivity-only when
used. Singleton groups make its group MCP penalty the ordinary coordinatewise
MCP penalty. The `grpreg` path does not expose a local-convexity index, so the
backward-compatible local-convexity field is `NA` and no such claim is made.

At thresholds above one half the package records the usual stability-selection
bound together with its assumptions. A threshold at or below one half is
labelled an exploratory frequency screen and receives no PFER claim.

Methodological foundations include Meinshausen and Buehlmann (2010,
doi:10.1111/j.1467-9868.2010.00740.x), Shah and Samworth (2013,
doi:10.1111/j.1467-9868.2011.01034.x), Belloni and Chernozhukov (2013,
doi:10.3150/11-BEJ410), and Dufour (2006,
doi:10.1016/j.jeconom.2005.06.007).

## Prespecified-core square-root-LASSO framework

Let

\[
y=X_A\beta_A+X_B\beta_B+\varepsilon,\qquad
H_0:\beta_B=0_q,
\]

where the prespecified core has \(p_1=|A|<n\) and the tested block may have
\(q=|B|\gg n\). Within each training fit, define

\[
X^{s}_{ij}=\frac{X_{ij}-\bar X_j}{s_j},\qquad
y_i^s=\frac{y_i-\bar y}{s_y},
\]

where \(s_j^2=n^{-1}\sum_i(X_{ij}-\bar X_j)^2\) and
\(s_y^2=n^{-1}\sum_i(y_i-\bar y)^2\). The primary full-model endpoint is one
joint all-predictor square-root LASSO:

\[
\widehat b^{FM}\in\arg\min_{b\in\mathbb R^p}
\left\{
\frac{\|y^s-X^sb\|_2}{\sqrt n}
+\lambda\sum_{j=1}^{p}\omega_j|b_j|
\right\},\qquad \omega_j>0.
\]

All \(p\) coordinates are estimated together and penalized; the core/tested
partition does not enter the FM optimization. The standardized intercept is
exactly zero. Returned coefficients and predictions are transformed back by
\[
\widehat\beta_j^{FM}=\frac{s_y}{s_j}\widehat b_j^{FM},\qquad
\widehat\alpha^{FM}=\bar y-\bar X^\top\widehat\beta^{FM}.
\]

The submodel endpoint is a fresh low-dimensional refit under the exact null:

\[
\widehat b^{SM}=((X_A^s)^+y^s,0_q).
\]

Thus SM uses all observations but only the prespecified \(X_A\) predictors and
sets the \(X_B\) coefficients to zero. It does not copy the core coordinates
from FM, and it is not supplied with the true active set. The core
pseudoinverse is computed from a checked economy SVD; no OLS Gram inverse
exists or is assumed. Its coefficients and intercept are back-transformed by
the same training-fit scaling rule.

For the Gaussian-iid design, the default test uses the exact-null residual.
Let \(M_A=I-X_A^s(X_A^s)^+\), \(r_A=M_Ay^s\),
\(d=n-1-\operatorname{rank}(X_A^s)\), and
\(z_j=M_AX^s_{B,j}\). Then

\[
T_{\mathrm{score}}=
\max_{j\in B}\frac{|z_j^\top r_A|}
{\|z_j\|_2\sqrt{\|r_A\|_2^2/d}}.
\]

The reported maximum partial-t statistic is the monotone transform

\[
T_{\max,t}=
\sqrt{\frac{(d-1)T_{\mathrm{score}}^2}
{d-T_{\mathrm{score}}^2}}.
\]

Gaussian draws are centered and projected off the same core, and the observed
and simulated maxima receive the same transformation. Thus the finite-draw
Monte Carlo p-value is invariant to the transformation. There is no direct
chi-square approximation, joint partial-F claim, or \(q\times q\) covariance
inverse. With

\[
\widehat\kappa=
\left\{\frac1B\sum_{b=1}^B(T_{\max,t}^{*(b)})^{-2}\right\}^{-1},
\]

the default estimators are

\[
\widehat\beta^{PT}
=\widehat\beta^{SM}
+\mathbf 1\{T_{\max,t}>c^*_{1-\alpha}\}
(\widehat\beta^{FM}-\widehat\beta^{SM}),
\]

\[
\widehat\beta^{S}
=\widehat\beta^{SM}
+\left(1-\frac{\widehat\kappa}
{T_{\max,t}^2\vee\varepsilon}\right)
(\widehat\beta^{FM}-\widehat\beta^{SM}),
\]

\[
\widehat\beta^{PS}
=\widehat\beta^{SM}
+\left[1-\frac{\widehat\kappa}
{T_{\max,t}^2\vee\varepsilon}\right]_+
(\widehat\beta^{FM}-\widehat\beta^{SM}).
\]

The max-calibrated Stein family is a research proposal. Classical
James--Stein dominance is not asserted. The legacy residual-multiplier test
and second-moment calibration remain available explicitly with
`test_calibration = "residual_multiplier"` and
`shrinkage_calibration = "second_moment"`.

## Minimal example

    set.seed(1)
    n <- 100
    p <- 1000
    p1 <- 20
    X <- matrix(rnorm(n * p), n, p)
    beta <- numeric(p)
    core_pattern <- c(
      1.50, 1.25, 1.00, 0.90, 0.80,
      -1.25, -1.00, -0.90, 0.75, -0.75
    )
    beta[1:20] <- rep(core_pattern, 2)
    y <- drop(X %*% beta + rnorm(n))

    fit <- fit_partial_sqrt_lasso_shrinkage(
      X,
      y,
      core_set = 1:p1,
      tested_set = (p1 + 1):p,
      full_endpoint = "all_x_square_root_lasso",
      all_x_penalty_loadings = rep(1, p),
      standardize_y = TRUE,
      bootstrap_B = 499,
      bootstrap_seed = 20260820
    )

    print(fit)
    coef(fit, method = "FM")
    coef(fit, method = "SM")
    coef(fit, method = "PS")

Here all 20 coordinates in the prespecified core are active. Under the exact
null, SM is the intended correct-restriction benchmark; the estimator is not
given the support separately and only receives the scientifically specified
core/tested partition.

The all-X FM endpoint does not require a diagonal population-precision
assertion. Conditional on fixed \(X\), the exact-null max-partial-\(t\) test
needs a prespecified core, nondegenerate residualized tested directions, and
homoskedastic Gaussian errors; it does not require independent columns of
\(X\).

## Risk convention

The accompanying simulations use the total coefficient loss

\[
L(\widehat\beta,\beta)=\|\widehat\beta-\beta\|_2^2
\]

and

\[
\operatorname{RPE}(m)=
\frac{\operatorname{MSE}(\widehat\beta^{FM})}
{\operatorname{MSE}(\widehat\beta^{m})}.
\]

With a fixed all-\(X\) Ridge penalty and full-vector coefficient loss, SM RPE
need not tend to zero.  Both Ridge FM and a misspecified SM can have
quadratic-in-signal loss because Ridge retains null-space bias when \(p>n\);
SM RPE then approaches a positive design-specific plateau.  The max-Stein
estimators should still approach FM and hence RPE one when their
data-dependent weight approaches one.  A sparse all-\(X\) LASSO sensitivity
may recover a strong sparse departure and produce an SM RPE close to zero,
but it is a different FM endpoint with additional assumptions.

## Retained research modes

- fit_ridge_mcp_shrinkage() retains the version 0.2.0 all-\(X\)
  Ridge--MCP sensitivity framework.
- fit_hd_shrinkage() retains the earlier square-root-LASSO restriction
  projection and small-\(q\) Wald branch for reproducibility.

The thresholded partial-debiased and profiled partial square-root-LASSO
endpoints are also retained as explicit backward-compatible modes. The
mandatory-core CPSS--SM/Ridge workflow is the version 0.6.0 data-adaptive
extension; the prespecified-core square-root-LASSO workflow and earlier modes remain
available for theory, sensitivity analysis, and reproducibility.
