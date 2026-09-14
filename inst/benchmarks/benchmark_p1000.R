#!/usr/bin/env Rscript

# Repository benchmark for the p = 1000, n = 100 Gaussian iid design.
# Run from the HD_shrinkage project root:
#   Rscript HDMaxShrink/inst/benchmarks/benchmark_p1000.R

package_dir <- normalizePath("HDMaxShrink", mustWork = TRUE)
prototype_file <- normalizePath(
  "normal_shrinkage_pggn.R", mustWork = TRUE
)
pkgload::load_all(package_dir, quiet = TRUE)
legacy <- new.env(parent = globalenv())
sys.source(prototype_file, envir = legacy)

set.seed(20260819)
n <- 100L
p <- 1000L
q <- 990L
B <- 499L
X_raw <- matrix(stats::rnorm(n * p), n, p)
y_raw <- stats::rnorm(n)
prep <- legacy$prepare_design(X_raw, y_raw)
lambda <- 1.1 * sqrt(2 * log(2 * p) / n)

time_r <- system.time({
  fit_r <- legacy$sqrt_lasso_primal_dual(
    prep$X,
    prep$y,
    lambda = lambda,
    max_iter = 30000L,
    tol = 5e-6,
    check_every = 50L
  )
})
time_cpp <- system.time({
  fit_cpp <- HDMaxShrink:::.sqrt_lasso_standardized(
    prep$X,
    prep$y,
    lambda = lambda,
    max_iter = 30000L,
    tol = 5e-6,
    check_every = 50L
  )
})

psi <- matrix(stats::rnorm(n * q), n, q)
psi <- sweep(psi, 2L, colMeans(psi), "-")
standard_error <- sqrt(colMeans(psi^2))
multipliers <- matrix(stats::rnorm(n * B), n, B)

time_boot_r <- system.time({
  bootstrap_r <- numeric(B)
  for (b in seq_len(B)) {
    draw <- as.vector(crossprod(multipliers[, b], psi)) /
      sqrt(n) / standard_error
    bootstrap_r[b] <- max(abs(draw))
  }
})
time_boot_cpp <- system.time({
  bootstrap_cpp <- HDMaxShrink:::cpp_multiplier_max(
    psi,
    standard_error,
    multipliers,
    block_size = 256L
  )
})

cat(sprintf(
  paste0(
    "n=%d p=%d q=%d B=%d\n",
    "solver_R=%.4fs solver_Cpp=%.4fs speedup=%.2fx\n",
    "solver_max_abs_diff=%.3e\n",
    "bootstrap_R=%.4fs bootstrap_Cpp=%.4fs speedup=%.2fx\n",
    "bootstrap_max_abs_diff=%.3e\n"
  ),
  n, p, q, B,
  time_r[["elapsed"]], time_cpp[["elapsed"]],
  time_r[["elapsed"]] / time_cpp[["elapsed"]],
  max(abs(fit_r$beta - fit_cpp$beta)),
  time_boot_r[["elapsed"]], time_boot_cpp[["elapsed"]],
  time_boot_r[["elapsed"]] / time_boot_cpp[["elapsed"]],
  max(abs(bootstrap_r - bootstrap_cpp))
))
