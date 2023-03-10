#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
double fit_ascr(NumericVector theta,
                arma::mat bin_capt,
                arma::mat mask,
                arma::mat traps,
                arma::mat mask_dists,
                arma::mat observed_bearings,
                arma::mat expected_bearings,
                double pixel_area) {
  // Transform the parameters back to their original scale
  double D = exp(theta[0]);
  double lambda0 = exp(theta[1]);
  double sigma = exp(theta[2]);
  double kappa = exp(theta[3]);

  // Set up the constants we will consistently work with
  arma::uword n_calls = bin_capt.n_rows;
  arma::uword n_traps = traps.n_rows;
  arma::uword n_mask = mask_dists.n_rows;

  // Evaluate the call detection function across the mask by trap
  arma::mat g_call(n_mask, n_traps);
  for (arma::uword j = 0; j < n_traps; j++) {
    g_call.col(j) = 1 - exp(-lambda0 * exp(-pow(mask_dists.col(j), 2) / (2 * pow(sigma, 2))));
  }

  // Evaluate the probability of detecting a call (at all) across the mask
  arma::rowvec p_call(n_mask);
  for (arma::uword s = 0; s < n_mask; s++) {
    p_call(s) = 1 - prod(1 - g_call.row(s)) + DBL_MIN;
  }

  // Compute the effective survey area
  double esa = pixel_area * sum(p_call);

  // Evaluate the log of the Poisson p.m.f modelling the number of detect calls
  double ell = R::dpois(n_calls, D * esa, 1);

  // Init f(s) & evaluate the filtered spatial point process
  arma::rowvec fs = (D * p_call) / (D * esa);

  // Init f(w|s)
  arma::mat fw(n_calls, n_mask);
  fw.ones();

  // Init f(y|w,s)
  arma::mat fy(n_calls, n_mask);
  fy.ones();

  // Evaluate the p.m.f./p.d.f. per pixel for f(w|s) and f(y|w,s)
  for (arma::uword i = 0; i < n_calls; i++) {
    for (arma::uword j = 0; j < n_traps; j++) {
      fw.row(i) %= pow(g_call.col(j).t(), bin_capt(i, j));
      fw.row(i) %= pow(1 - g_call.col(j).t(), 1 - bin_capt(i, j));

      fy.row(i) %= pow(exp(kappa * cos(observed_bearings(i, j) - expected_bearings.col(j).t())), bin_capt(i, j));
      fy.row(i) /= pow(2 * M_PI * R::bessel_i(kappa, 0, 1), bin_capt(i, j));
    }
    fw.row(i) /= p_call;

    // Then add each call's contribution to the log-likelihood
    ell += log(pixel_area * sum(fy.row(i) % fw.row(i) % fs + DBL_MIN));
  }

  Rcout << "nll: " << -ell << " D: " << D << " lambda0: " << lambda0 << " sigma: " << sigma << " kappa: " << kappa << " esa: " << esa << " \n";

  // Return the negative log-likelihood
  return -ell;
}
