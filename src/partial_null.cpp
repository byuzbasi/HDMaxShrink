#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

//' SVD geometry for a prespecified low-dimensional core model
//'
//' The routine never forms `(X_core' X_core)^{-1}`.  It returns the
//' Moore--Penrose coefficient operator (after verifying full column rank), an
//' orthonormal basis for the core column space, and every tested column
//' residualized against that space.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_partial_null_geometry(const arma::mat& X_core,
                                     const arma::mat& X_tested,
                                     const double rank_tolerance = 1e-10) {
  const arma::uword n = X_core.n_rows;
  const arma::uword p1 = X_core.n_cols;
  const arma::uword q = X_tested.n_cols;

  if (n == 0 || p1 == 0 || p1 >= n || !X_core.is_finite()) {
    Rcpp::stop(
      "X_core must be finite, non-empty, and have fewer columns than rows."
    );
  }
  if (X_tested.n_rows != n || q == 0 || !X_tested.is_finite()) {
    Rcpp::stop(
      "X_tested must be finite, non-empty, and have nrow(X_core) rows."
    );
  }
  if (!std::isfinite(rank_tolerance) || rank_tolerance <= 0.0) {
    Rcpp::stop("rank_tolerance must be finite and positive.");
  }

  arma::mat U;
  arma::vec singular_values;
  arma::mat V;
  const bool ok = arma::svd_econ(U, singular_values, V, X_core);
  if (!ok || singular_values.n_elem != p1 ||
      !singular_values.is_finite() || singular_values[0] <= 0.0) {
    Rcpp::stop("The economy SVD of X_core failed.");
  }

  const double threshold = rank_tolerance *
    static_cast<double>(std::max(n, p1)) * singular_values[0];
  const arma::uword numerical_rank = arma::accu(singular_values > threshold);
  if (numerical_rank != p1) {
    Rcpp::stop(
      "X_core is not full column rank at the requested tolerance."
    );
  }

  arma::mat scaled_V = V;
  scaled_V.each_row() /= singular_values.t();
  const arma::mat coefficient_operator = scaled_V * U.t();
  arma::mat residualized_tested = X_tested - U * (U.t() * X_tested);
  const arma::vec residualized_second_moment =
    arma::sum(arma::square(residualized_tested), 0).t() /
    static_cast<double>(n);
  if (!residualized_second_moment.is_finite() ||
      arma::any(residualized_second_moment <= 0.0)) {
    Rcpp::stop("Every residualized tested column must have positive variance.");
  }
  const double orthogonality_error =
    arma::abs(U.t() * residualized_tested).max();
  const double condition_number =
    singular_values[0] / singular_values[singular_values.n_elem - 1];

  return Rcpp::List::create(
    Rcpp::Named("coefficient_operator") = coefficient_operator,
    Rcpp::Named("orthonormal_basis") = U,
    Rcpp::Named("residualized_tested") = residualized_tested,
    Rcpp::Named("residualized_second_moment") =
      residualized_second_moment,
    Rcpp::Named("singular_values") = singular_values,
    Rcpp::Named("rank") = static_cast<int>(numerical_rank),
    Rcpp::Named("rank_tolerance") = rank_tolerance,
    Rcpp::Named("condition_number") = condition_number,
    Rcpp::Named("orthogonality_error") = orthogonality_error
  );
}


//' One-step scores for the profiled partial square-root LASSO
//'
//' The tested design has already been residualized against the prespecified
//' core.  A diagonal population precision is used, which is appropriate for
//' the Gaussian-iid design studied in the accompanying simulations.  No
//' tested-block covariance matrix is formed or inverted.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_partial_debiased_scores(
    const arma::vec& y_residualized,
    const arma::mat& residualized_tested,
    const arma::vec& beta_tested) {
  const arma::uword n = y_residualized.n_elem;
  const arma::uword q = residualized_tested.n_cols;

  if (n == 0 || q == 0 || !y_residualized.is_finite() ||
      residualized_tested.n_rows != n ||
      !residualized_tested.is_finite() ||
      beta_tested.n_elem != q || !beta_tested.is_finite()) {
    Rcpp::stop("The partial one-step inputs have incompatible dimensions.");
  }

  const arma::vec second_moment =
    arma::sum(arma::square(residualized_tested), 0).t() /
    static_cast<double>(n);
  if (!second_moment.is_finite() || arma::any(second_moment <= 0.0)) {
    Rcpp::stop("Every residualized tested column must have positive variance.");
  }

  const arma::vec residual =
    y_residualized - residualized_tested * beta_tested;
  const arma::vec correction =
    (residualized_tested.t() * residual) /
    static_cast<double>(n) / second_moment;
  const arma::vec theta_tilde = beta_tested + correction;

  arma::mat psi = residualized_tested;
  psi.each_col() %= residual;
  psi.each_row() /= second_moment.t();
  const arma::rowvec psi_mean = arma::mean(psi, 0);
  psi.each_row() -= psi_mean;
  const arma::vec variance =
    arma::sum(arma::square(psi), 0).t() / static_cast<double>(n);
  if (!variance.is_finite() || arma::any(variance <= 0.0)) {
    Rcpp::stop("Every partial one-step score must have positive variance.");
  }

  return Rcpp::List::create(
    Rcpp::Named("theta_tilde") = theta_tilde,
    Rcpp::Named("correction") = correction,
    Rcpp::Named("residual") = residual,
    Rcpp::Named("psi_centered") = psi,
    Rcpp::Named("variance") = variance,
    Rcpp::Named("second_moment") = second_moment,
    Rcpp::Named("residual_norm") = arma::norm(residual, 2)
  );
}


//' Apply cached core-model geometry to a centered response
//'
//' The returned influence scores are the heteroskedasticity-robust,
//' residualized score contributions used by the multiplier-bootstrap maximum
//' test of the coordinate null.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_partial_null_apply(const arma::vec& y,
                                  const arma::mat& coefficient_operator,
                                  const arma::mat& orthonormal_basis,
                                  const arma::mat& residualized_tested) {
  const arma::uword n = y.n_elem;
  const arma::uword p1 = coefficient_operator.n_rows;
  const arma::uword q = residualized_tested.n_cols;

  if (n == 0 || p1 == 0 || q == 0 || !y.is_finite()) {
    Rcpp::stop("The partial-null inputs must be finite and non-empty.");
  }
  if (coefficient_operator.n_cols != n ||
      orthonormal_basis.n_rows != n ||
      orthonormal_basis.n_cols != p1 ||
      residualized_tested.n_rows != n ||
      !coefficient_operator.is_finite() ||
      !orthonormal_basis.is_finite() ||
      !residualized_tested.is_finite()) {
    Rcpp::stop("The cached partial-null geometry has incompatible dimensions.");
  }

  const arma::vec beta_core = coefficient_operator * y;
  const arma::vec residual =
    y - orthonormal_basis * (orthonormal_basis.t() * y);
  arma::mat psi = residualized_tested;
  psi.each_col() %= residual;
  const arma::rowvec score_mean = arma::mean(psi, 0);
  psi.each_row() -= score_mean;
  const arma::vec variance =
    arma::sum(arma::square(psi), 0).t() / static_cast<double>(n);
  if (!variance.is_finite() || arma::any(variance <= 0.0)) {
    Rcpp::stop("Every residualized score must have positive finite variance.");
  }

  return Rcpp::List::create(
    Rcpp::Named("beta_core") = beta_core,
    Rcpp::Named("residual") = residual,
    Rcpp::Named("score_mean") = score_mean.t(),
    Rcpp::Named("psi_centered") = psi,
    Rcpp::Named("variance") = variance,
    Rcpp::Named("residual_norm") = arma::norm(residual, 2)
  );
}
