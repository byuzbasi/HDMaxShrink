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
  old <- par(mar = c(4, 4.6, 3.8, 1.2), family = "sans", las = 1,
             fg = "#39434D", col.axis = "#39434D", col.lab = "#39434D")
  on.exit(par(old))
  upper <- max(results$Test_MSE) * 1.20
  positions <- barplot(
    results$Test_MSE, names.arg = results$Method,
    col = c("#84939D", "#009E73", "#E69F00", "#56B4E9", "#0072B2"),
    border = NA, ylim = c(0, upper), axes = FALSE,
    ylab = "Independent test MSE", cex.names = 1.1
  )
  axis(2, lwd = 0, lwd.ticks = 0, cex.axis = 0.95)
  text(positions, results$Test_MSE, labels = sprintf("%.3f", results$Test_MSE),
       pos = 3, cex = 0.95)
  title(main = "One held-out sample, five estimators", adj = 0, cex.main = 1.15)
  mtext("p = 1,000 | selection n = 100 | analysis n = 100 | test n = 200",
        side = 3, line = 0.5, adj = 0, cex = 0.82, col = "#65717C")
  invisible(results)
}
if (interactive()) plot_quickstart()
cat("Small example passed; one draw is not evidence of uniform dominance.\n")
