# Small reproducible example only; no production simulation.
library(HDMaxShrink)

set.seed(20260914)
n <- 60L
p <- 120L
beta <- c(1.5, -1.25, 1, rep(0, p - 3L))
make_sample <- function() {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("x", seq_len(p))
  list(X = X, y = drop(X %*% beta + rnorm(n)))
}
selection_data <- make_sample()
analysis_data <- make_sample()

# Required predictors are fixed before examining either sample.
selection <- cpss_select_core(
  selection_data$X, selection_data$y,
  selector = "lasso", mandatory_core = 1:3,
  complementary_pairs = 3L, base_selection_size = 3L,
  stability_threshold = 0.60, path_points = 30L, seed = 20260915L
)
fit <- fit_cpss_ridge_shrinkage(
  analysis_data$X, analysis_data$y,
  selection = selection, ridge_lambda = 0.25,
  bootstrap_B = 99L, bootstrap_seed = 20260916L
)
methods <- c("FM", "SM", "PT", "S", "PS")
# coef() includes the back-transformed intercept as its first entry.
estimates <- vapply(methods, function(m) coef(fit, method = m), numeric(p + 1L))
stopifnot(identical(dim(estimates), c(p + 1L, 5L)), all(is.finite(estimates)))
print(round(estimates[1:6, , drop = FALSE], 4))
cat("Small example passed; this is not a performance comparison.\n")
