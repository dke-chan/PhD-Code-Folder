#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
double fit_ascrStack(NumericVector theta,
                     arma::mat bin_capt,
                     arma::mat mask,
                     arma::mat traps,
                     arma::mat mask_dists,
                     arma::mat observed_bearings,
                     arma::mat expected_bearings,
                     double pixel_area,
                     bool bernoulli,
                     arma::vec ID,
                     double survey_length) {
  // Transform the parameters back to their original scale
  double D = exp(theta[0]);
  double lambda0 = exp(theta[1]);
  double sigma = exp(theta[2]);
  double callProcessTheta = theta[3];
  if (bernoulli) {
    callProcessTheta = R::plogis(callProcessTheta, 0, 1, 1, 0);
  } else {
    callProcessTheta = exp(callProcessTheta);
  }
  double kappa = exp(theta[4]);

  // Set up the constants we will consistently work with
  arma::uword n_calls = bin_capt.n_rows;
  arma::uword n_traps = traps.n_rows;
  arma::uword n_mask = mask_dists.n_rows;
  arma::vec animal_ID = unique(ID);
  arma::uword n_animal = animal_ID.n_elem;
  arma::vec n_indiv_calls(n_animal);
  n_indiv_calls.zeros();

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

  // Construct the elements for f(c|s)
  arma::mat fc_lmnt(2, n_mask);
  arma::rowvec p_animal(n_mask);

  if (bernoulli) {
    fc_lmnt.row(0) = 1 - callProcessTheta * p_call;
    fc_lmnt.row(1) = callProcessTheta * p_call;

    // Evaluate the probability of detecting an animal (at all ~ assuming stacking) across the mask
    p_animal = 1 - pow(fc_lmnt.row(0), survey_length) + DBL_MIN;
  } else {
    fc_lmnt.row(0) = 1 - exp(-callProcessTheta * p_call);
    fc_lmnt.row(1) = callProcessTheta * p_call;

    // Evaluate the probability of detecting an animal (at all ~ assuming stacking) across the mask
    p_animal = fc_lmnt.row(0) + DBL_MIN;
  }

  // Compute the effective survey area
  double esa = pixel_area * sum(p_animal);

  // Evaluate the log of the Poisson p.m.f modelling the number of detect calls
  double ell = R::dpois(n_animal, D * esa, 1);

  // Init f(s) & evaluate the filtered spatial point process
  arma::rowvec fs = (D * p_animal) / (D * esa);

  // Init f(O|c,s)
  arma::mat fO(n_animal, n_mask);
  fO.ones();

  // Init f(Y|O,c,s)
  arma::mat fY(n_animal, n_mask);
  fY.ones();

  // Evaluate the p.m.f./p.d.f. per pixel for f(O|c,s) and f(Y|O,c,s)
  for (arma::uword j = 0; j < n_calls; j++) {
    arma::uword i = ID(j);
    for (arma::uword k = 0; k < n_traps; k++) {
      fO.row(i) %= pow(g_call.col(k).t(), bin_capt(j, k));
      fO.row(i) %= pow(1 - g_call.col(k).t(), 1 - bin_capt(j, k));

      fY.row(i) %= pow(exp(kappa * cos(observed_bearings(j, k) - expected_bearings.col(k).t())), bin_capt(j, k));
      fY.row(i) /= pow(2 * M_PI * R::bessel_i(kappa, 0, 1), bin_capt(j, k));
    }
    fO.row(i) /= p_call;

    // Also increment the number of calls that animal made by one
    n_indiv_calls(i) += 1;
  }

  for (arma::uword i = 0; i < n_animal; i++) {
    // Init f(c|s) & evaluate the p.m.f. per pixel for f(c|s)
    arma::rowvec fc(n_mask);
    if (bernoulli) {
      //fc = pow(fc_lmnt.row(1), n_indiv_calls(i)) % pow(fc_lmnt.row(0), survey_length - n_indiv_calls(i)) / p_animal;
      fc = R::choose(survey_length, n_indiv_calls(i)) * pow(fc_lmnt.row(1), n_indiv_calls(i)) % pow(fc_lmnt.row(0), survey_length - n_indiv_calls(i)) / p_animal;
    } else {
      fc = (pow(fc_lmnt.row(1), n_indiv_calls(i)) % exp(-callProcessTheta * p_call) / tgamma(n_indiv_calls(i) + 1)) % (1 / p_animal);
    }

    // Then add each animal's contribution to the log-likelihood
    ell += log(pixel_area * sum(fY.row(i) % fO.row(i) % fc % fs + DBL_MIN));
  }

  if (bernoulli) {
    Rcout << "nll: " << -ell << " D: " << D << " lambda0: " << lambda0 << " sigma: " << sigma << " pBernoulli: " << callProcessTheta << " kappa: " << kappa << " esa: " << esa << "\n";
  } else {
    Rcout << "nll: " << -ell << " D: " << D << " lambda0: " << lambda0 << " sigma: " << sigma << " lambdaC: " << callProcessTheta << " kappa: " << kappa << " esa: " << esa << "\n";
  }

  return -ell;
}

// [[Rcpp::export]]
double fit_ascrStackSession(NumericVector theta, arma::field<arma::mat> bin_capt, arma::field<arma::mat> mask, arma::field<arma::mat> traps,
                            arma::field<arma::mat> mask_dists, arma::field<arma::mat> observed_bearings, arma::field<arma::mat> expected_bearings,
                            arma::vec pixel_area, bool bernoulli, arma::field<arma::vec> ID, arma::vec survey_length) {
  // Transform the parameters back to their original scale
  double D = exp(theta[0]);
  double lambda0 = exp(theta[1]);
  double sigma = exp(theta[2]);
  double callProcessTheta = theta[3];
  if (bernoulli) {
    callProcessTheta = R::plogis(callProcessTheta, 0, 1, 1, 0);
  } else {
    callProcessTheta = exp(callProcessTheta);
  }
  double kappa = exp(theta[4]);

  // Set up the constants we will consistently work with
  arma::uword n_sessions = bin_capt.n_elem;
  arma::vec n_calls(n_sessions);
  arma::vec n_traps(n_sessions);
  arma::vec n_mask(n_sessions);
  arma::vec n_animal(n_sessions);

  for (arma::uword e = 0; e < n_sessions; e++) {
    n_calls(e) = bin_capt(e).n_rows;
    n_traps(e) = traps(e).n_rows;
    n_mask(e) = mask(e).n_rows;
    arma::vec animal_ID = unique(ID(e));
    n_animal(e) = animal_ID.n_elem;
  }

  // ...
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

    // Construct the elements for f(c|s)
    arma::mat fc_lmnt(n_mask(e), 2);
    arma::vec p_animal(n_mask(e));

    if (bernoulli) {
      fc_lmnt.col(0) = 1 - callProcessTheta * p_call;
      fc_lmnt.col(1) = callProcessTheta * p_call;

      // Evaluate the probability of detecting an animal (at all ~ assuming stacking) across the mask
      p_animal = 1 - pow(fc_lmnt.col(0), survey_length(e)) + DBL_MIN;
    } else {
      fc_lmnt.col(0) = 1 - exp(-callProcessTheta * p_call);
      fc_lmnt.col(1) = callProcessTheta * p_call;

      // Evaluate the probability of detecting an animal (at all ~ assuming stacking) across the mask
      p_animal = fc_lmnt.col(0) + DBL_MIN;
    }

    // Compute the effective survey area
    double esa_e = pixel_area(e) * sum(p_animal);
    esa += esa_e;

    // Init f(s) & evaluate the filtered spatial point process
    arma::vec fs = (D * p_animal) / (D * esa_e);

    // Init f(O|c,s)
    arma::mat fO(n_animal(e), n_mask(e));
    fO.ones();

    // Init f(Y|O,c,s)
    arma::mat fY(n_animal(e), n_mask(e));
    fY.ones();

    // ...
    arma::vec n_indiv_calls(n_animal(e));
    n_indiv_calls.zeros();

    // Evaluate the p.m.f./p.d.f. per pixel for f(O|c,s) and f(Y|O,c,s)
    for (arma::uword j = 0; j < n_calls(e); j++) {
      arma::uword i = ID(e)(j);
      for (arma::uword k = 0; k < n_traps(e); k++) {
        fO.row(i) %= pow(g_call.col(k).t(), bin_capt(e)(j, k));
        fO.row(i) %= pow(1 - g_call.col(k).t(), 1 - bin_capt(e)(j, k));

        fY.row(i) %= pow(exp(kappa * cos(observed_bearings(e)(j, k) - expected_bearings(e).col(k).t())), bin_capt(e)(j, k));
        fY.row(i) /= pow(2 * M_PI * R::bessel_i(kappa, 0, 1), bin_capt(e)(j, k));
      }
      fO.row(i) /= p_call.t();

      // Also increment the number of calls that animal made by one
      n_indiv_calls(i) += 1;
    }

    for (arma::uword i = 0; i < n_animal(e); i++) {
      // Init f(c|s) & evaluate the p.m.f. per pixel for f(c|s)
      arma::vec fc(n_mask(e));
      if (bernoulli) {
        fc = R::choose(survey_length(e), n_indiv_calls(i)) * pow(fc_lmnt.col(1), n_indiv_calls(i)) % pow(fc_lmnt.col(0), survey_length(e) - n_indiv_calls(i)) / p_animal;
      } else {
        fc = (pow(fc_lmnt.col(1), n_indiv_calls(i)) % exp(-callProcessTheta * p_call) / tgamma(n_indiv_calls(i) + 1)) % (1 / p_animal);
      }

      // Then add each animal's contribution to the log-likelihood
      ell += log(pixel_area(e) * sum(fY.row(i) % fO.row(i) % fc.t() % fs.t()) + DBL_MIN);
    }
  }

  // Evaluate the log of the Poisson p.m.f modelling the number of detect calls
  ell += R::dpois(sum(n_animal), D * esa, 1);

  if (bernoulli) {
    Rcout << "nll: " << -ell << " D: " << D << " lambda0: " << lambda0 << " sigma: " << sigma << " pBernoulli: " << callProcessTheta << " kappa: " << kappa << " esa: " << esa << "\n";
  } else {
    Rcout << "nll: " << -ell << " D: " << D << " lambda0: " << lambda0 << " sigma: " << sigma << " lambdaC: " << callProcessTheta << " kappa: " << kappa << " esa: " << esa << "\n";
  }

  return -ell;
}