#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

namespace hdms {

double combined_operator_norm(const arma::mat& X,
                              const arma::mat& A,
                              const int max_iter,
                              const double tol) {
  const arma::uword p = X.n_cols;
  arma::vec z(p);
  for (arma::uword j = 0; j < p; ++j) {
    const double index = static_cast<double>(j + 1);
    z[j] = std::sin(index) + std::cos(index / 3.0);
  }
  const double initial_norm = arma::norm(z, 2);
  if (!std::isfinite(initial_norm) || initial_norm == 0.0) {
    Rcpp::stop("Could not initialize the operator-norm iteration.");
  }
  z /= initial_norm;

  double old_value = 0.0;
  double value = 0.0;
  for (int iter = 0; iter < max_iter; ++iter) {
    const arma::vec Xz = X * z;
    arma::vec direction = X.t() * Xz;
    double squared_norm = arma::dot(Xz, Xz);

    if (A.n_rows > 0) {
      const arma::vec Az = A * z;
      direction += A.t() * Az;
      squared_norm += arma::dot(Az, Az);
    }

    const double direction_norm = arma::norm(direction, 2);
    if (direction_norm == 0.0) {
      return 0.0;
    }
    z = direction / direction_norm;
    value = squared_norm;

    if (std::abs(value - old_value) <=
        tol * std::max(1.0, value)) {
      break;
    }
    old_value = value;
  }

  // Match the safety inflation used by the R research prototype.
  return 1.02 * std::sqrt(value);
}

arma::vec soft_threshold(const arma::vec& value,
                         const arma::vec& threshold) {
  arma::vec result(value.n_elem);
  for (arma::uword j = 0; j < value.n_elem; ++j) {
    const double magnitude = std::abs(value[j]) - threshold[j];
    if (magnitude <= 0.0) {
      result[j] = 0.0;
    } else {
      result[j] = std::copysign(magnitude, value[j]);
    }
  }
  return result;
}

double kkt_residual_l1(const arma::vec& beta,
                       const arma::vec& gradient,
                       const double lambda,
                       const arma::vec& penalty_factor,
                       const double zero_tol = 1e-8) {
  double residual = 0.0;
  for (arma::uword j = 0; j < beta.n_elem; ++j) {
    double current = 0.0;
    if (std::abs(beta[j]) > zero_tol) {
      current = std::abs(
        gradient[j] + lambda * penalty_factor[j] *
        std::copysign(1.0, beta[j])
      );
    } else {
      current = std::max(
        std::abs(gradient[j]) - lambda * penalty_factor[j], 0.0
      );
    }
    residual = std::max(residual, current);
  }
  return residual;
}

}  // namespace hdms


//' RcppArmadillo operator norm for the stacked data/restriction operator
//'
//' @keywords internal
// [[Rcpp::export]]
double cpp_operator_norm(const arma::mat& X,
                         const arma::mat& A,
                         const int max_iter = 200,
                         const double tol = 1e-9) {
  if (X.n_rows == 0 || X.n_cols == 0 || !X.is_finite()) {
    Rcpp::stop("X must be a non-empty finite matrix.");
  }
  if (A.n_cols != X.n_cols || !A.is_finite()) {
    Rcpp::stop("A must be finite and have ncol(X) columns.");
  }
  if (max_iter < 1 || !std::isfinite(tol) || tol <= 0.0) {
    Rcpp::stop("Invalid operator-norm iteration controls.");
  }
  return hdms::combined_operator_norm(X, A, max_iter, tol);
}


//' RcppArmadillo primal-dual square-root lasso solver
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_sqrt_lasso_pd(const arma::mat& X,
                             const arma::vec& y,
                             const double lambda,
                             const arma::mat& A,
                             const arma::vec& target,
                             const arma::vec& penalty_factor,
                             const arma::vec& beta_init,
                             const int max_iter = 30000,
                             const double tol = 5e-6,
                             const int check_every = 50,
                             const bool verbose = false,
                             const double operator_norm_value = NA_REAL) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  const arma::uword q = A.n_rows;

  if (n == 0 || p == 0 || !X.is_finite()) {
    Rcpp::stop("X must be a non-empty finite matrix.");
  }
  if (y.n_elem != n || !y.is_finite()) {
    Rcpp::stop("y must be finite and have nrow(X) entries.");
  }
  if (!std::isfinite(lambda) || lambda <= 0.0) {
    Rcpp::stop("lambda must be finite and positive.");
  }
  if (A.n_cols != p || !A.is_finite()) {
    Rcpp::stop("A must be finite and have ncol(X) columns.");
  }
  if (target.n_elem != q || !target.is_finite()) {
    Rcpp::stop("target must be finite and have nrow(A) entries.");
  }
  if (penalty_factor.n_elem != p || !penalty_factor.is_finite() ||
      arma::any(penalty_factor < 0.0)) {
    Rcpp::stop("penalty_factor must contain p finite nonnegative values.");
  }
  if (beta_init.n_elem != p || !beta_init.is_finite()) {
    Rcpp::stop("beta_init must contain p finite values.");
  }
  if (max_iter < 1 || check_every < 1 || !std::isfinite(tol) || tol <= 0.0) {
    Rcpp::stop("Invalid solver iteration controls.");
  }

  double L = operator_norm_value;
  if (Rcpp::NumericVector::is_na(L)) {
    L = hdms::combined_operator_norm(X, A, 200, 1e-9);
  }
  if (!std::isfinite(L) || L <= 0.0) {
    Rcpp::stop("The combined design/restriction operator has zero norm.");
  }

  const double primal_step = 0.95 / L;
  const double dual_step = 0.95 / L;
  const double dual_radius = 1.0 / std::sqrt(static_cast<double>(n));

  arma::vec beta = beta_init;
  arma::vec beta_bar = beta;
  arma::vec dual_data(n, arma::fill::zeros);
  arma::vec dual_restriction(q, arma::fill::zeros);

  bool converged = false;
  double objective = std::numeric_limits<double>::infinity();
  double relative_change = std::numeric_limits<double>::infinity();
  double restriction_violation = q > 0 ?
    std::numeric_limits<double>::infinity() : 0.0;
  double kkt_residual = std::numeric_limits<double>::infinity();
  int completed_iterations = 0;

  for (int iter = 1; iter <= max_iter; ++iter) {
    completed_iterations = iter;
    const arma::vec beta_old = beta;

    dual_data += dual_step * (X * beta_bar - y);
    const double dual_norm = arma::norm(dual_data, 2);
    if (!std::isfinite(dual_norm)) {
      Rcpp::stop("Non-finite dual norm encountered.");
    }
    if (dual_norm > dual_radius && dual_norm > 0.0) {
      dual_data *= dual_radius / dual_norm;
    }

    if (q > 0) {
      dual_restriction += dual_step * (A * beta_bar - target);
    }

    arma::vec gradient = X.t() * dual_data;
    if (q > 0) {
      gradient += A.t() * dual_restriction;
    }
    const arma::vec threshold =
      primal_step * lambda * penalty_factor;
    beta = hdms::soft_threshold(
      beta_old - primal_step * gradient, threshold
    );
    beta_bar = 2.0 * beta - beta_old;

    if (iter % check_every == 0 || iter == max_iter) {
      const arma::vec residual = y - X * beta;
      objective = arma::norm(residual, 2) /
        std::sqrt(static_cast<double>(n)) +
        lambda * arma::dot(penalty_factor, arma::abs(beta));
      relative_change = arma::norm(beta - beta_old, 2) /
        std::max(1.0, arma::norm(beta_old, 2));

      arma::vec stationarity_gradient = X.t() * dual_data;
      if (q > 0) {
        stationarity_gradient += A.t() * dual_restriction;
        restriction_violation = arma::abs(A * beta - target).max();
      } else {
        restriction_violation = 0.0;
      }
      kkt_residual = hdms::kkt_residual_l1(
        beta, stationarity_gradient, lambda, penalty_factor
      );

      double target_scale = 1.0;
      if (q > 0) {
        target_scale = std::max(1.0, arma::abs(target).max());
      }
      const double scaled_restriction =
        restriction_violation / target_scale;
      if (std::max({relative_change, kkt_residual, scaled_restriction}) <=
          tol) {
        converged = true;
        break;
      }

      if (verbose && iter % (20 * check_every) == 0) {
        Rcpp::Rcout << "iter=" << iter
                    << " objective=" << objective
                    << " rel=" << relative_change
                    << " kkt=" << kkt_residual
                    << " restriction=" << restriction_violation
                    << "\n";
      }
    }

    if (iter % 256 == 0) {
      Rcpp::checkUserInterrupt();
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("iterations") = completed_iterations,
    Rcpp::Named("objective") = objective,
    Rcpp::Named("relative_change") = relative_change,
    Rcpp::Named("kkt_residual") = kkt_residual,
    Rcpp::Named("restriction_violation") = restriction_violation,
    Rcpp::Named("operator_norm") = L,
    Rcpp::Named("lambda") = lambda
  );
}
