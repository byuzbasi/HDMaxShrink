#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

namespace hdms {

Rcpp::List finalize_scores(const arma::vec& theta_tilde,
                           const arma::vec& residual,
                           arma::mat psi) {
  const arma::uword n = psi.n_rows;
  if (n == 0 || psi.n_cols == 0) {
    Rcpp::stop("The influence-score matrix must be non-empty.");
  }
  const arma::rowvec means = arma::mean(psi, 0);
  psi.each_row() -= means;
  const arma::vec variance =
    arma::sum(arma::square(psi), 0).t() / static_cast<double>(n);
  if (!variance.is_finite() || arma::any(variance <= 0.0)) {
    Rcpp::stop("Every target must have a finite positive score variance.");
  }
  return Rcpp::List::create(
    Rcpp::Named("theta_tilde") = theta_tilde,
    Rcpp::Named("residual") = residual,
    Rcpp::Named("psi_centered") = psi,
    Rcpp::Named("variance") = variance
  );
}

}  // namespace hdms


//' Debiased target scores for general precision directions
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_debiased_scores_dense(const arma::mat& X,
                                     const arma::vec& y,
                                     const arma::vec& beta_full,
                                     const arma::vec& target_point,
                                     const arma::mat& directions) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  const arma::uword q = directions.n_rows;
  if (n == 0 || p == 0 || !X.is_finite()) {
    Rcpp::stop("X must be a non-empty finite matrix.");
  }
  if (y.n_elem != n || beta_full.n_elem != p ||
      target_point.n_elem != q || directions.n_cols != p ||
      !y.is_finite() || !beta_full.is_finite() ||
      !target_point.is_finite() || !directions.is_finite()) {
    Rcpp::stop("Incompatible or non-finite debiasing inputs.");
  }

  const arma::vec residual = y - X * beta_full;
  const arma::vec empirical_gradient =
    X.t() * residual / static_cast<double>(n);
  const arma::vec theta_tilde =
    target_point + directions * empirical_gradient;
  arma::mat psi = X * directions.t();
  psi.each_col() %= residual;
  return hdms::finalize_scores(theta_tilde, residual, std::move(psi));
}


//' Debiased coordinate scores for a diagonal precision model
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_debiased_scores_diagonal(const arma::mat& X,
                                        const arma::vec& y,
                                        const arma::vec& beta_full,
                                        const Rcpp::IntegerVector& targets,
                                        const arma::vec& direction_scale,
                                        const arma::vec& target_point) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  const arma::uword q = targets.size();
  if (n == 0 || p == 0 || !X.is_finite()) {
    Rcpp::stop("X must be a non-empty finite matrix.");
  }
  if (y.n_elem != n || beta_full.n_elem != p ||
      direction_scale.n_elem != q || target_point.n_elem != q ||
      !y.is_finite() || !beta_full.is_finite() ||
      !direction_scale.is_finite() || !target_point.is_finite()) {
    Rcpp::stop("Incompatible or non-finite diagonal debiasing inputs.");
  }

  arma::uvec index(q);
  for (arma::uword k = 0; k < q; ++k) {
    const int current = targets[k];
    if (current < 1 || current > static_cast<int>(p)) {
      Rcpp::stop("targets must use one-based indices in 1,...,p.");
    }
    index[k] = static_cast<arma::uword>(current - 1);
  }

  const arma::vec residual = y - X * beta_full;
  const arma::vec target_gradient =
    X.cols(index).t() * residual / static_cast<double>(n);
  const arma::vec theta_tilde =
    target_point + direction_scale % target_gradient;
  arma::mat psi = X.cols(index);
  psi.each_row() %= direction_scale.t();
  psi.each_col() %= residual;
  return hdms::finalize_scores(theta_tilde, residual, std::move(psi));
}


//' Multiplier-bootstrap maxima from centered influence scores
//'
//' @keywords internal
// [[Rcpp::export]]
arma::vec cpp_multiplier_max(const arma::mat& psi_centered,
                             const arma::vec& standard_error,
                             const arma::mat& multipliers,
                             const int block_size = 256) {
  const arma::uword n = psi_centered.n_rows;
  const arma::uword q = psi_centered.n_cols;
  const arma::uword draws = multipliers.n_cols;
  if (n == 0 || q == 0 || draws == 0 ||
      multipliers.n_rows != n ||
      standard_error.n_elem != q ||
      !psi_centered.is_finite() || !multipliers.is_finite() ||
      !standard_error.is_finite() || arma::any(standard_error <= 0.0)) {
    Rcpp::stop("Invalid multiplier-bootstrap inputs.");
  }
  if (block_size < 1) {
    Rcpp::stop("block_size must be positive.");
  }

  arma::vec result(draws);
  const double root_n = std::sqrt(static_cast<double>(n));
  for (arma::uword first = 0; first < draws;
       first += static_cast<arma::uword>(block_size)) {
    const arma::uword last = std::min(
      draws - 1,
      first + static_cast<arma::uword>(block_size) - 1
    );
    arma::mat standardized =
      psi_centered.t() * multipliers.cols(first, last) / root_n;
    standardized.each_col() /= standard_error;
    const arma::rowvec block_max = arma::max(arma::abs(standardized), 0);
    result.subvec(first, last) = block_max.t();
    Rcpp::checkUserInterrupt();
  }
  return result;
}


//' Conditional Gaussian score maxima for an exact coordinate null
//'
//' The input residual must be the centered exact-null residual, the tested
//' design must already be residualized off the centered core, and the columns
//' of `core_basis` must be an orthonormal basis for that core. Each Gaussian
//' draw is centered and projected off the core before applying the same
//' residual-scale studentization as the observed statistic.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_gaussian_score_max(const arma::vec& null_residual,
                                  const arma::mat& residualized_tested,
                                  const arma::mat& core_basis,
                                  const arma::mat& gaussian_draws,
                                  const int block_size = 256) {
  const arma::uword n = residualized_tested.n_rows;
  const arma::uword q = residualized_tested.n_cols;
  const arma::uword rank = core_basis.n_cols;
  const arma::uword draws = gaussian_draws.n_cols;
  if (n < 2 || q == 0 || draws == 0 ||
      null_residual.n_elem != n || core_basis.n_rows != n ||
      gaussian_draws.n_rows != n ||
      !null_residual.is_finite() || !residualized_tested.is_finite() ||
      !core_basis.is_finite() || !gaussian_draws.is_finite()) {
    Rcpp::stop("Invalid conditional Gaussian score inputs.");
  }
  if (rank + 1 >= n) {
    Rcpp::stop("The centered null residual degrees of freedom must be positive.");
  }
  if (block_size < 1) {
    Rcpp::stop("block_size must be positive.");
  }

  const arma::vec direction_norm = arma::sqrt(
    arma::sum(arma::square(residualized_tested), 0).t()
  );
  if (!direction_norm.is_finite() || arma::any(direction_norm <= 0.0)) {
    Rcpp::stop("Every residualized tested column must have positive norm.");
  }
  const double residual_df = static_cast<double>(n - rank - 1);
  const double sigma_hat = std::sqrt(
    arma::dot(null_residual, null_residual) / residual_df
  );
  if (!std::isfinite(sigma_hat) || sigma_hat <= 0.0) {
    Rcpp::stop("The exact-null residual scale must be finite and positive.");
  }

  const arma::vec observed_scores =
    residualized_tested.t() * null_residual / direction_norm / sigma_hat;
  const double observed = arma::abs(observed_scores).max();
  arma::vec bootstrap_statistics(draws);

  for (arma::uword first = 0; first < draws;
       first += static_cast<arma::uword>(block_size)) {
    const arma::uword last = std::min(
      draws - 1,
      first + static_cast<arma::uword>(block_size) - 1
    );
    arma::mat residual_draws = gaussian_draws.cols(first, last);
    const arma::rowvec means = arma::mean(residual_draws, 0);
    residual_draws.each_row() -= means;
    if (rank > 0) {
      residual_draws -= core_basis * (core_basis.t() * residual_draws);
    }
    const arma::rowvec draw_scale = arma::sqrt(
      arma::sum(arma::square(residual_draws), 0) / residual_df
    );
    if (!draw_scale.is_finite() || arma::any(draw_scale <= 0.0)) {
      Rcpp::stop("A Gaussian null draw has invalid residual scale.");
    }
    arma::mat standardized = residualized_tested.t() * residual_draws;
    standardized.each_col() /= direction_norm;
    standardized.each_row() /= draw_scale;
    bootstrap_statistics.subvec(first, last) =
      arma::max(arma::abs(standardized), 0).t();
    Rcpp::checkUserInterrupt();
  }

  return Rcpp::List::create(
    Rcpp::Named("statistic") = observed,
    Rcpp::Named("T_max") = observed,
    Rcpp::Named("score") = observed_scores,
    Rcpp::Named("sigma_hat") = sigma_hat,
    Rcpp::Named("residual_df") = residual_df,
    Rcpp::Named("direction_norm") = direction_norm,
    Rcpp::Named("bootstrap_statistics") = bootstrap_statistics
  );
}
