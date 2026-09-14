#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

namespace hdms {

double signed_value(const double magnitude, const double reference) {
  return reference >= 0.0 ? magnitude : -magnitude;
}

double mcp_coordinate_update(const double z,
                             const double curvature,
                             const double lambda,
                             const double gamma) {
  const double absolute_z = std::abs(z);
  if (absolute_z <= lambda) {
    return 0.0;
  }
  if (absolute_z <= gamma * lambda * curvature) {
    const double denominator = curvature - 1.0 / gamma;
    if (denominator <= 0.0) {
      Rcpp::stop(
        "MCP coordinate curvature must exceed the inverse concavity parameter."
      );
    }
    return signed_value((absolute_z - lambda) / denominator, z);
  }
  return z / curvature;
}

double mcp_penalty(const double coefficient,
                   const double lambda,
                   const double gamma) {
  const double magnitude = std::abs(coefficient);
  if (magnitude <= gamma * lambda) {
    return lambda * magnitude - magnitude * magnitude / (2.0 * gamma);
  }
  return 0.5 * gamma * lambda * lambda;
}

}  // namespace hdms


//' Dual high-dimensional ridge path
//'
//' Computes ridge estimates from an `n` by `n` eigendecomposition rather than
//' a `p` by `p` inverse. The objective is
//' `||y-X beta||^2/(2n) + lambda ||beta||^2/2`.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ridge_dual_path(const arma::mat& X,
                               const arma::vec& y,
                               const arma::vec& lambda) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  const arma::uword path_length = lambda.n_elem;
  if (n == 0 || p == 0 || !X.is_finite()) {
    Rcpp::stop("X must be a non-empty finite matrix.");
  }
  if (y.n_elem != n || !y.is_finite()) {
    Rcpp::stop("y must be finite and have nrow(X) entries.");
  }
  if (path_length == 0 || !lambda.is_finite() || arma::any(lambda <= 0.0)) {
    Rcpp::stop("lambda must contain finite positive values.");
  }

  const arma::mat kernel = X * X.t();
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  if (!arma::eig_sym(eigenvalues, eigenvectors, kernel)) {
    Rcpp::stop("The ridge kernel eigendecomposition failed.");
  }
  eigenvalues.transform([](double value) {
    return std::max(0.0, value);
  });
  const arma::vec rotated_y = eigenvectors.t() * y;

  arma::mat beta(p, path_length, arma::fill::zeros);
  arma::vec rss(path_length, arma::fill::zeros);
  arma::vec degrees_freedom(path_length, arma::fill::zeros);
  arma::vec gcv(path_length, arma::fill::zeros);

  for (arma::uword index = 0; index < path_length; ++index) {
    const arma::vec denominator = eigenvalues +
      static_cast<double>(n) * lambda[index];
    const arma::vec alpha = eigenvectors * (rotated_y / denominator);
    beta.col(index) = X.t() * alpha;
    const arma::vec residual = y - X * beta.col(index);
    rss[index] = arma::dot(residual, residual);
    degrees_freedom[index] = arma::sum(eigenvalues / denominator);
    // The original-scale model includes one unpenalized intercept. Although
    // X and y are centered before this kernel, GCV must count that intercept
    // in the total effective degrees of freedom.
    const double total_degrees_freedom = degrees_freedom[index] + 1.0;
    const double leverage_gap = std::max(
      1e-12,
      1.0 - total_degrees_freedom / static_cast<double>(n)
    );
    gcv[index] = (rss[index] / static_cast<double>(n)) /
      (leverage_gap * leverage_gap);
  }

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta,
    Rcpp::Named("lambda") = lambda,
    Rcpp::Named("rss") = rss,
    Rcpp::Named("degrees_freedom") = degrees_freedom,
    Rcpp::Named("gcv") = gcv,
    Rcpp::Named("kernel_eigenvalues") = eigenvalues,
    Rcpp::Named("converged") = true
  );
}


//' MCP coordinate-descent path on all design columns
//'
//' The objective is `||y-X beta||^2/(2n)` plus the coordinatewise MCP
//' penalty. Successive decreasing penalties use warm starts.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_mcp_path(const arma::mat& X,
                        const arma::vec& y,
                        const arma::vec& lambda,
                        const arma::vec& beta_init,
                        const double gamma = 3.0,
                        const int max_iter = 5000,
                        const double tol = 1e-7,
                        const double zero_tol = 1e-8) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  const arma::uword path_length = lambda.n_elem;
  if (n == 0 || p == 0 || !X.is_finite()) {
    Rcpp::stop("X must be a non-empty finite matrix.");
  }
  if (y.n_elem != n || !y.is_finite()) {
    Rcpp::stop("y must be finite and have nrow(X) entries.");
  }
  if (path_length == 0 || !lambda.is_finite() || arma::any(lambda <= 0.0)) {
    Rcpp::stop("lambda must contain finite positive values.");
  }
  for (arma::uword index = 1; index < path_length; ++index) {
    if (lambda[index] > lambda[index - 1]) {
      Rcpp::stop("The MCP lambda path must be nonincreasing.");
    }
  }
  if (!std::isfinite(gamma) || gamma <= 1.0) {
    Rcpp::stop("gamma must be finite and exceed one.");
  }
  if (max_iter < 1 || !std::isfinite(tol) || tol <= 0.0 ||
      !std::isfinite(zero_tol) || zero_tol < 0.0) {
    Rcpp::stop("Invalid MCP iteration controls.");
  }

  arma::vec coefficient;
  if (beta_init.n_elem == 0) {
    coefficient.zeros(p);
  } else {
    if (beta_init.n_elem != p || !beta_init.is_finite()) {
      Rcpp::stop("beta_init must be empty or contain p finite values.");
    }
    coefficient = beta_init;
  }
  const arma::vec curvature =
    arma::sum(arma::square(X), 0).t() / static_cast<double>(n);
  if (!curvature.is_finite() || arma::any(curvature <= 1.0 / gamma)) {
    Rcpp::stop(
      "Every standardized column curvature must exceed 1/gamma for MCP."
    );
  }

  arma::mat beta_path(p, path_length, arma::fill::zeros);
  arma::vec rss(path_length, arma::fill::zeros);
  arma::vec objective(path_length, arma::fill::zeros);
  arma::ivec degrees_freedom(path_length, arma::fill::zeros);
  Rcpp::LogicalVector converged(path_length);
  Rcpp::IntegerVector iterations(path_length);
  arma::vec maximum_change(path_length, arma::fill::zeros);
  arma::vec residual = y - X * coefficient;

  for (arma::uword path_index = 0;
       path_index < path_length;
       ++path_index) {
    bool path_converged = false;
    double last_change = std::numeric_limits<double>::infinity();
    int completed_iterations = 0;

    for (int iteration = 1; iteration <= max_iter; ++iteration) {
      completed_iterations = iteration;
      double max_change = 0.0;
      for (arma::uword column = 0; column < p; ++column) {
        const double old_value = coefficient[column];
        if (old_value != 0.0) {
          residual += X.col(column) * old_value;
        }
        const double z = arma::dot(X.col(column), residual) /
          static_cast<double>(n);
        double new_value = hdms::mcp_coordinate_update(
          z, curvature[column], lambda[path_index], gamma
        );
        if (std::abs(new_value) <= zero_tol) {
          new_value = 0.0;
        }
        coefficient[column] = new_value;
        if (new_value != 0.0) {
          residual -= X.col(column) * new_value;
        }
        max_change = std::max(
          max_change,
          std::abs(new_value - old_value) * std::sqrt(curvature[column])
        );
      }
      last_change = max_change;
      if (max_change <= tol * std::max(1.0, arma::norm(coefficient, 2))) {
        path_converged = true;
        break;
      }
      if (iteration % 64 == 0) {
        Rcpp::checkUserInterrupt();
      }
    }

    beta_path.col(path_index) = coefficient;
    rss[path_index] = arma::dot(residual, residual);
    double penalty = 0.0;
    for (arma::uword column = 0; column < p; ++column) {
      penalty += hdms::mcp_penalty(
        coefficient[column], lambda[path_index], gamma
      );
    }
    objective[path_index] = rss[path_index] /
      (2.0 * static_cast<double>(n)) + penalty;
    degrees_freedom[path_index] = static_cast<int>(
      arma::accu(arma::abs(coefficient) > zero_tol)
    );
    converged[path_index] = path_converged;
    iterations[path_index] = completed_iterations;
    maximum_change[path_index] = last_change;
    Rcpp::checkUserInterrupt();
  }

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta_path,
    Rcpp::Named("lambda") = lambda,
    Rcpp::Named("gamma") = gamma,
    Rcpp::Named("rss") = rss,
    Rcpp::Named("objective") = objective,
    Rcpp::Named("degrees_freedom") = degrees_freedom,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("maximum_change") = maximum_change
  );
}
