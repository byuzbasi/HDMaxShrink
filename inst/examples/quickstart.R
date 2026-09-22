# One reproducible toy dataset; no production simulation or risk estimate.
library(HDMaxShrink)

# 1. Generate three independent samples from the same sparse Gaussian model.
set.seed(20260914)
n <- 100L
p <- 1000L
n_test <- 200L
beta <- c(1.50, 1.25, 1.00, 0.90, 0.80,
          -1.25, -1.00, -0.90, 0.75, -0.75, rep(0, p - 10L))
make_sample <- function(rows = n) {
  X <- matrix(rnorm(rows * p), nrow = rows, ncol = p)
  colnames(X) <- paste0("x", seq_len(p))
  list(X = X, y = drop(X %*% beta + rnorm(rows)))
}
selection_data <- make_sample()
analysis_data <- make_sample()
test_data <- make_sample(n_test)

# 2. Only x1 and x2 are assumed known in advance; CPSS sees no true beta.
# Remaining signals x3:x10 must be learned from the selection sample.
# Twenty complementary pairs give 40 half-sample fits (50 rows each).
# The budget caps each half-fit, NOT the size of the final submodel.
selection <- cpss_select_core(
  selection_data$X, selection_data$y,
  selector = "lasso", mandatory_core = 1:2,
  complementary_pairs = 20L, base_selection_size = 10L,
  stability_threshold = 0.60, path_points = 30L, seed = 20260915L
)

# 3. Fit on independent analysis data. Lambda is fixed, not test-set tuned.
# The package centers/scales analysis X and y and back-transforms predictions.
fit <- fit_cpss_ridge_shrinkage(
  analysis_data$X, analysis_data$y,
  selection = selection, ridge_lambda = 0.25,
  bootstrap_B = 199L, bootstrap_seed = 20260916L,
  shrinkage_calibration = "inverse_moment"
)
methods <- c("FM", "SM", "PT", "S", "PS")
# coef() includes the back-transformed intercept as its first entry.
estimates <- vapply(methods, function(m) coef(fit, method = m), numeric(p + 1L))

# 4. Report what was selected and how the test defines the interpolation.
core_table <- selection$stability_table[
  selection$stability_table$selected_for_core,
  c("feature", "core_role", "stability_frequency"), drop = FALSE
]
cat("Selected submodel (NA frequency means prespecified, not selected):\n")
print(core_table, row.names = FALSE)
# Mandatory variables have no CPSS frequency: they are included by design.
optional_table <- selection$stability_table[
  is.finite(selection$stability_table$stability_frequency), , drop = FALSE
]
optional_table$half_sample_selections <- as.integer(round(
  optional_table$stability_frequency * 2L * selection$complementary_pairs
))
cat("Top optional candidates (ranked on selection data only):\n")
print(head(optional_table[, c("feature", "half_sample_selections",
                             "stability_frequency", "selected_for_core")], 10L),
      row.names = FALSE)

# A nonzero intercept is not a selected predictor. All vectors return p slopes.
coefficient_summary <- data.frame(
  Method = methods,
  Slopes_returned = p,
  Nonzero_slopes = colSums(abs(estimates[-1L, , drop = FALSE]) > 1e-8),
  Intercept = unname(estimates[1L, ]), row.names = NULL
)
cat("Coefficient counts (numerical tolerance 1e-8, not significance tests):\n")
print(coefficient_summary, row.names = FALSE, digits = 5)
test_summary <- data.frame(
  n = n, p = p, p1 = length(fit$core_set), q = length(fit$tested_set),
  T_max = fit$inference$statistic, p_value = fit$inference$p_value,
  reject = fit$reject, kappa = fit$inference$shrinkage_calibration
)
print(test_summary, row.names = FALSE, digits = 5)

# Weights multiply (FM - SM): 0 is SM, 1 is FM. S may be negative.
full_weight <- c(FM = 1, SM = 0, PT = as.numeric(fit$reject),
                 S = fit$shrinkage$stein_weight,
                 PS = fit$shrinkage$positive_weight)

# 5. Evaluate once on untouched test data, with the training scaling map.
predictions <- vapply(methods, function(m) {
  predict(fit, newdata = test_data$X, method = m)
}, numeric(n_test))
results_table <- data.frame(
  Method = methods,
  Full_weight = unname(full_weight[methods]),
  Coefficient_loss = colSums((estimates[-1L, , drop = FALSE] - beta)^2),
  Test_MSE = colMeans((predictions - test_data$y)^2),
  row.names = NULL
)
cat("Single-dataset results (coefficient loss excludes the intercept):\n")
print(results_table, row.names = FALSE, digits = 5)
stopifnot(
  identical(dim(estimates), c(p + 1L, 5L)),
  identical(dim(predictions), c(n_test, 5L)),
  all(is.finite(estimates)), all(is.finite(predictions)),
  all(is.finite(as.matrix(results_table[, -1L]))),
  all(1:2 %in% fit$core_set),
  all(estimates[-1L, "SM"][fit$tested_set] == 0)
)

# 6. Plot on the current graphics device; no files are written automatically.
plot_quickstart <- function(results = results_table) {
  old <- par(no.readonly = TRUE)
  on.exit({ layout(1); par(old) })
  colors <- c(FM = "#8793A0", SM = "#007F6E", PT = "#C48212",
              S = "#64A9CB", PS = "#215BB5")
  ink <- "#243347"
  muted <- "#5B6B7C"
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1.25, 1))
  par(oma = c(2.7, 0.6, 4.8, 0.6), family = "sans", las = 1,
      fg = ink, col.axis = muted, col.lab = ink, cex = 0.92)

  # Display ranking uses selection data, not truth or held-out performance.
  shown <- head(optional_table, 10L)
  position <- rev(seq_len(nrow(shown)))
  selected <- shown$selected_for_core
  par(mar = c(3.8, 4.6, 4.1, 2.1))
  plot(NA, xlim = c(0, 0.83), ylim = c(0.5, nrow(shown) + 0.5),
       axes = FALSE, xlab = "Selection frequency across 40 half-samples", ylab = "",
       xaxs = "i", yaxs = "i")
  abline(v = seq(0, 0.8, 0.2), col = "#E8EDF2", lwd = 0.8)
  segments(0, position, shown$stability_frequency, position,
           col = ifelse(selected, "#91CFC1", "#CED7E1"), lwd = 3)
  abline(v = selection$stability_threshold, col = colors["SM"], lty = 2, lwd = 1.2)
  points(shown$stability_frequency, position, pch = 21, cex = 1.15,
         bg = ifelse(selected, colors["SM"], "#AFBBC7"), col = "white", lwd = 1)
  axis(1, at = seq(0, 0.8, 0.2), labels = sprintf("%.1f", seq(0, 0.8, 0.2)),
       lwd = 0, lwd.ticks = 0)
  axis(2, at = position, labels = shown$feature, lwd = 0, lwd.ticks = 0)
  text(rep(0.71, nrow(shown)), position,
       labels = paste0(shown$half_sample_selections, "/40"), adj = 0,
       cex = 0.78, col = ifelse(selected, colors["SM"], muted))
  title("A  |  CPSS learns the extension", adj = 0, cex.main = 1.06, line = 2.6)
  mtext("Top 10 of 998 optional predictors; threshold = 0.60", side = 3,
        line = 1.4, adj = 0, cex = 0.74, col = muted)
  mtext("x1 + x2: mandatory, outside the frequency screen", side = 3,
        line = 0.3, adj = 0, cex = 0.74, col = colors["SM"])

  par(mar = c(3.8, 4.0, 4.1, 2.2))
  position <- rev(seq_len(nrow(results)))
  upper <- max(results$Test_MSE) * 1.22
  plot(NA, xlim = c(0, upper), ylim = c(0.5, nrow(results) + 0.5),
       axes = FALSE, xlab = "Independent test MSE  (lower is better)", ylab = "",
       xaxs = "i", yaxs = "i")
  ticks <- pretty(c(0, upper), n = 5)
  ticks <- ticks[ticks >= 0 & ticks <= upper]
  abline(v = ticks, col = "#E8EDF2", lwd = 0.8)
  segments(0, position, results$Test_MSE, position,
           col = unname(colors[results$Method]), lwd = 5)
  points(results$Test_MSE, position, pch = 21, cex = 1.5,
         bg = unname(colors[results$Method]), col = "white", lwd = 1.2)
  axis(1, at = ticks, lwd = 0, lwd.ticks = 0)
  axis(2, at = position, labels = results$Method, lwd = 0, lwd.ticks = 0)
  text(results$Test_MSE + upper * 0.045, position,
       sprintf("%.3f", results$Test_MSE), adj = 0, cex = 0.88, col = ink)
  title("B  |  Evaluate on untouched data", adj = 0, cex.main = 1.06, line = 2.6)
  mtext("200 test observations; same rows for every method", side = 3,
        line = 1.4, adj = 0, cex = 0.74, col = muted)
  mtext("Ridge FM: 1,000 predictors   |   Selected SM: 3", side = 3,
        line = 0.3, adj = 0, cex = 0.74, col = muted)

  # PT=FM and S=PS in this draw; avoid hiding overplotted duplicate curves.
  shown_indices <- seq_len(20L)
  par(mar = c(4.7, 4.6, 3.9, 2.2))
  plot(NA, xlim = c(0.6, 20.4), ylim = c(-1.65, 1.85), axes = FALSE,
       xlab = "Predictor index (first 20 shown; all 1,000 enter the losses)",
       ylab = "Coefficient", xaxs = "i")
  rect(10.5, -1.65, 20.4, 1.85, col = "#F3F6FA", border = NA)
  abline(h = c(-1, 0, 1), col = c("#E5EBF1", "#B3BFCC", "#E5EBF1"), lwd = 0.8)
  for (m in c("FM", "SM", "PS")) {
    lines(shown_indices, estimates[shown_indices + 1L, m],
          col = colors[m], lwd = 1.3, type = "b", pch = c(FM = 1, SM = 15, PS = 17)[m],
          cex = 0.68)
  }
  points(shown_indices, beta[shown_indices], pch = 4, col = ink, cex = 1, lwd = 1.5)
  axis(1, at = shown_indices, lwd = 0, lwd.ticks = 0, cex.axis = 0.86)
  axis(2, at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5), lwd = 0, lwd.ticks = 0)
  title("C  |  See what selection and shrinkage retain", adj = 0, cex.main = 1.06, line = 2.5)
  legend("topright", legend = c("Truth", "FM = PT", "SM", "S = PS"),
         col = c(ink, colors[c("FM", "SM", "PS")]), pch = c(4, 1, 15, 17),
         lty = c(NA, 1, 1, 1), lwd = c(1.5, 1.3, 1.3, 1.3), bty = "n",
         ncol = 4, cex = 0.82, x.intersp = 1.4,
         text.width = strwidth("FM = PT", cex = 0.82) * 1.25,
         inset = c(0, -0.23), xpd = NA)
  mtext("CPSS selection  >  Ridge / submodel  >  Max-test shrinkage", outer = TRUE,
        side = 3, line = 2.7, adj = 0.035, cex = 1.25, font = 2, col = ink)
  mtext("Selection n = 100  |  Analysis n = 100  |  Test n = 200  |  p = 1,000",
        outer = TRUE, side = 3, line = 1.1, adj = 0.035, cex = 0.88, col = muted)
  mtext("Fixed-seed illustration, not a simulation study. No guarantee of support recovery or risk dominance.",
        outer = TRUE, side = 1, line = 1.1, adj = 0.035, cex = 0.77, col = muted)
  invisible(results)
}
if (interactive()) plot_quickstart()
cat("Small example passed; one draw is not evidence of uniform dominance.\n")
