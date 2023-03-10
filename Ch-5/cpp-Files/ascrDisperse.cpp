#include <RcppArmadillo.h>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::vec dbinorm(arma::mat x, arma::rowvec mu, double sigma) {
  arma::vec xout = (x.col(0) - mu(0)) / sigma;
  arma::vec yout = (x.col(1) - mu(1)) / sigma;

  return (exp(-1 * (pow(xout, 2) + pow(yout, 2)) / 2) / (2 * M_PI * sigma * sigma));
}
// https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Bivariate_case

arma::uvec whichElems(arma::mat s, arma::rowvec t, double constraint) {
  arma::vec euclideanDists = sqrt(pow(s.col(0) - t(0), 2) + pow(s.col(1) - t(1), 2));

  return (find(euclideanDists < constraint));
}
//

// [[Rcpp::export]]
double fit_ascrDisperse(NumericVector theta,
                        arma::mat bin_capt,
                        arma::mat mask,
                        arma::mat traps,
                        arma::mat mask_dists,
                        arma::mat observed_bearings,
                        arma::mat expected_bearings,
                        double pixel_area,
                        arma::vec ID,
                        double survey_length,
                        double threshold = 1e-4,
                        bool trace = false) {
  // Transform the parameters back to their original scale
  double D = exp(theta[0]);
  double lambda0 = exp(theta[1]);
  double sigma = exp(theta[2]);
  double pBinomial = R::plogis(theta[3], 0, 1, 1, 0);
  // double pBinomial = exp(theta[3]); // Actually lambdaC
  double sigmaMove = exp(theta[4]);
  double kappa = exp(theta[5]);

  // Set up the constants we will consistently work with
  arma::uword n_calls = bin_capt.n_rows;
  arma::uword n_traps = traps.n_rows;
  arma::uword n_mask = mask_dists.n_rows;
  arma::vec animal_ID = unique(ID);
  arma::uword n_animal = animal_ID.n_elem;
  arma::vec n_indiv_calls(n_animal);
  n_indiv_calls.zeros();
  double m2_pixel_area = pixel_area * 10000;
  double m2_pixel_edge = sqrt(m2_pixel_area);

  // Evaluate the call detection function across the mask by trap
  arma::mat g_call(n_mask, n_traps);
  for (arma::uword j = 0; j < n_traps; j++) {
    g_call.col(j) = (1 - exp(-lambda0 * exp(-pow(mask_dists.col(j), 2) / (2 * pow(sigma, 2)))));
  }

  // Evaluate the probability of detecting a call at location s (at all) across the mask
  arma::vec p_call(n_mask);
  for (arma::uword s = 0; s < n_mask; s++) {
    p_call(s) = 1 - prod(1 - g_call.row(s)) + DBL_MIN;
  }

  arma::mat fsStarDense(n_mask, n_mask);
  arma::vec p_actCen(n_mask);

  double biNormConstraint = R::qnorm5((2 - threshold) / 2, 0, 1, 1, 0) * sigmaMove;

  for (arma::uword t = 0; t < n_mask; t++) {
    arma::uvec indices = whichElems(mask, mask.row(t), biNormConstraint);

    arma::vec fs(n_mask);
    fs.zeros();
    if (biNormConstraint > m2_pixel_edge) {
      fs(indices) = m2_pixel_area * p_call(indices) % dbinorm(mask.rows(indices), mask.row(t), sigmaMove);
    } else {
      // If the biNormContraint is less than m2_pixel_edge, this probably means that the whole pixel encapsulates the binormal pdf...
      fs(indices) = p_call(indices);
    }

    // Evaluate the numerator for the filtered Bivariate Normal per pixel for f(s|t)
    fsStarDense.col(t) = fs;

    // Evaluate the probability of detecting a call at location s given location t (at all) across the mask
    // This is also the denominator for the filtered Bivariate Normal per pixel for f(s|t)
    p_actCen(t) = sum(fsStarDense.col(t));

    // Compute the filtered Bivariate Normal
    fsStarDense.col(t) /= p_actCen(t);
  }

  arma::sp_mat fsStar(fsStarDense);

  // Evaluate the probability of detecting an animal (at all) across the mask
  arma::vec p_animal = 1 - pow(1 - pBinomial * p_actCen, survey_length) + DBL_MIN; // Binomial
  // arma::vec p_animal = 1 - exp(-pBinomial * p_call_t) + DBL_MIN; // Poisson

  // Compute the effective survey area
  double esa = pixel_area * sum(p_animal);

  // Evaluate the log of the Poisson p.m.f modelling the number of detect calls
  double ell = R::dpois(n_animal, D * esa, 1);

  // Init f(t) & evaluate the filtered spatial point process
  arma::vec ft = (D * p_animal) / (D * esa);

  // Init f(Y,O|c,t)
  arma::mat fYO_t(n_animal, n_mask);
  fYO_t.ones();

  // Transpose p_call
  arma::rowvec p_call_transpose = p_call.t();

  // Evaluate the p.m.f./p.d.f. per pixel for f(O|c,s) and f(Y|O,c,s)
  for (arma::uword j = 0; j < n_calls; j++) {
    arma::uword i = ID(j);
    arma::rowvec fyos(n_mask);
    fyos.ones();

    for (arma::uword k = 0; k < n_traps; k++) {
      fyos %= pow(g_call.col(k).t(), bin_capt(j, k));
      fyos %= pow(1 - g_call.col(k).t(), 1 - bin_capt(j, k));

      fyos %= pow(exp(kappa * cos(observed_bearings(j, k) - expected_bearings.col(k).t())), bin_capt(j, k));
      fyos /= pow(2 * M_PI * R::bessel_i(kappa, 0, 1), bin_capt(j, k));
    }
    fyos /= p_call_transpose;

    // Also increment the number of calls that animal made by one
    n_indiv_calls(i) += 1;

    // Update the animal's contribution to the likelihood assuming independence between detections
    fYO_t.row(i) %= fyos * fsStar;
  }

  for (arma::uword i = 0; i < n_animal; i++) {
    // Evaluate the p.m.f. per pixel for f(c|s)
    arma::vec fc = R::choose(survey_length, n_indiv_calls(i)) * pow(pBinomial * p_actCen, n_indiv_calls(i)) % pow(1 - pBinomial * p_actCen, survey_length - n_indiv_calls(i)) / p_animal;
    // Binomial
    // arma::vec fc = (pow(pBinomial * p_call_t, n_indiv_calls(i)) % exp(-pBinomial * p_call_t) / R::gammafn(n_indiv_calls(i) + 1)) / p_animal;
    // Poisson

    // Then add each animal's contribution to the log-likelihood
    ell += log(sum(fYO_t.row(i).t() % fc % ft) + DBL_MIN);
  }

  if (trace) {
    Rcout << "nll: " << -ell << " D: " << D << " lambda0: " << lambda0 << " sigma: " << sigma << " pBinomial: " << pBinomial << " sigmaMove: " << sigmaMove << " kappa: " << kappa << " esa_na: " << esa << " esa_nc: " << pixel_area * sum(p_call) << "\n";
  }

  return -ell;
}

// [[Rcpp::export]]
double fit_ascrDisperseSession(NumericVector theta, arma::field<arma::mat> bin_capt,
                               arma::field<arma::mat> mask, arma::field<arma::mat> traps,
                               arma::field<arma::mat> mask_dists, arma::field<arma::mat> observed_bearings,
                               arma::field<arma::mat> expected_bearings, arma::vec pixel_area,
                               arma::field<arma::vec> ID, arma::vec survey_length,
                               double threshold = 1e-4, bool trace = false) {
  // Transform the parameters back to their original scale
  double D = exp(theta[0]);
  double lambda0 = exp(theta[1]);
  double sigma = exp(theta[2]);
  double pBinomial = R::plogis(theta[3], 0, 1, 1, 0);
  // double pBinomial = exp(theta[3]); // Actually lambdaC
  double sigmaMove = exp(theta[4]);
  double kappa = exp(theta[5]);

  // Set up the constants we will consistently work with
  arma::uword n_sessions = bin_capt.n_elem;
  arma::vec n_calls(n_sessions);
  arma::vec n_traps(n_sessions);
  arma::vec n_mask(n_sessions);
  arma::vec n_animal(n_sessions);
  arma::vec m2_pixel_area = pixel_area * 10000;
  arma::vec m2_pixel_edge = sqrt(m2_pixel_area);

  for (arma::uword e = 0; e < n_sessions; e++) {
    n_calls(e) = bin_capt(e).n_rows;
    n_traps(e) = traps(e).n_rows;
    n_mask(e) = mask(e).n_rows;
    arma::vec animal_ID = unique(ID(e));
    n_animal(e) = animal_ID.n_elem;
  }

  // arma::vec n_indiv_calls(n_animal);
  // n_indiv_calls.zeros();

  // ...
  double biNormConstraint = R::qnorm5((2 - threshold) / 2, 0, 1, 1, 0) * sigmaMove;
  double ell = 0;
  double esa = 0;

  // For each session...
  for (arma::uword e = 0; e < n_sessions; e++) {
    // Evaluate the call detection function across the mask by trap
    arma::mat g_call(n_mask(e), n_traps(e));
    for (arma::uword j = 0; j < n_traps(e); j++) {
      g_call.col(j) = (1 - exp(-lambda0 * exp(-pow(mask_dists(e).col(j), 2) / (2 * pow(sigma, 2)))));
    }

    // Evaluate the probability of detecting a call at location s (at all) across the mask
    arma::vec p_call(n_mask(e));
    for (arma::uword s = 0; s < n_mask(e); s++) {
      p_call(s) = 1 - prod(1 - g_call.row(s)) + DBL_MIN;
    }

    // ...
    arma::mat fsStarDense(n_mask(e), n_mask(e));
    arma::vec p_actCen(n_mask(e));

    for (arma::uword t = 0; t < n_mask(e); t++) {
      arma::uvec indices = whichElems(mask(e), mask(e).row(t), biNormConstraint);

      arma::vec fs(n_mask(e));
      fs.zeros();
      if (biNormConstraint > m2_pixel_edge(e)) {
        fs(indices) = m2_pixel_area(e) * p_call(indices) % dbinorm(mask(e).rows(indices), mask(e).row(t), sigmaMove);
      } else {
        fs(indices) = p_call(indices);
      }

      // Evaluate the numerator for the filtered Bivariate Normal per pixel for f(s|t)
      fsStarDense.col(t) = fs;

      // Evaluate the probability of detecting a call at location s given location t (at all) across the mask
        // This is also the denominator for the filtered Bivariate Normal per pixel for f(s|t)
      p_actCen(t) = sum(fsStarDense.col(t));

      // Compute the filtered Bivariate Normal
      fsStarDense.col(t) /= p_actCen(t);
    }

    arma::sp_mat fsStar(fsStarDense);

    // Evaluate the probability of detecting an animal (at all) across the mask
    arma::vec p_animal = 1 - pow(1 - pBinomial * p_actCen, survey_length(e)) + DBL_MIN; // Binomial
    // arma::vec p_animal = 1 - exp(-pBinomial * p_call_t) + DBL_MIN; // Poisson

    // Compute the effective survey area
    double esa_e = pixel_area(e) * sum(p_animal);
    esa += esa_e;

    // Init f(t) & evaluate the filtered spatial point process
    arma::vec ft = (D * p_animal) / (D * esa_e);

    // Init f(Y,O|c,t)
    arma::mat fYO_t(n_animal(e), n_mask(e));
    fYO_t.ones();

    // ...
    arma::vec n_indiv_calls(n_animal(e));
    n_indiv_calls.zeros();

    // Transpose p_call
    arma::rowvec p_call_transpose = p_call.t();

    // Evaluate the p.m.f./p.d.f. per pixel for f(O|c,s) and f(Y|O,c,s)
    for (arma::uword j = 0; j < n_calls(e); j++) {
      arma::uword i = ID(e)(j);
      arma::rowvec fyos(n_mask(e));
      fyos.ones();

      for (arma::uword k = 0; k < n_traps(e); k++) {
        fyos %= pow(g_call.col(k).t(), bin_capt(e)(j, k));
        fyos %= pow(1 - g_call.col(k).t(), 1 - bin_capt(e)(j, k));

        fyos %= pow(exp(kappa * cos(observed_bearings(e)(j, k) - expected_bearings(e).col(k).t())), bin_capt(e)(j, k));
        fyos /= pow(2 * M_PI * R::bessel_i(kappa, 0, 1), bin_capt(e)(j, k));
      }
      fyos /= p_call_transpose;

      // Also increment the number of calls that animal made by one
      n_indiv_calls(i) += 1;

      // Update the animal's contribution to the likelihood assuming independence between detections
      fYO_t.row(i) %= fyos * fsStar;
    }

    for (arma::uword i = 0; i < n_animal(e); i++) {
      // Evaluate the p.m.f. per pixel for f(c|s)
      arma::vec fc = R::choose(survey_length(e), n_indiv_calls(i)) * pow(pBinomial * p_actCen, n_indiv_calls(i)) % pow(1 - pBinomial * p_actCen, survey_length(e) - n_indiv_calls(i)) / p_animal;
      // Binomial
      // arma::vec fc = (pow(pBinomial * p_call_t, n_indiv_calls(i)) % exp(-pBinomial * p_call_t) / R::gammafn(n_indiv_calls(i) + 1)) / p_animal;
      // Poisson

      // Then add each animal's contribution to the log-likelihood
      ell += log(sum(fYO_t.row(i).t() % fc % ft) + DBL_MIN);
    }
  }

  // Evaluate the log of the Poisson p.m.f modelling the number of detect calls
  ell += R::dpois(sum(n_animal), D * esa, 1);

  if (trace) {
    Rcout << "nll: " << -ell << " D: " << D << " lambda0: " << lambda0 << " sigma: " << sigma << " pBinomial: " << pBinomial << " sigmaMove: " << sigmaMove << " kappa: " << kappa << "\n";
  }

  return -ell;
}