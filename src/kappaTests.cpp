/*** Resampling Cohen's kappa tests */

#include <Rcpp.h>
#include <Rmath.h>
#include <random>

#include "numericUtils.h"
#include "vectorUtils.h"
#include "transformUtils.h"
#include "rankUtils.h"
#include "covariance.h"
#include "kappa.h"

using namespace Rcpp;

// permutation test for a pair of vectors

//[[Rcpp::export]]

NumericVector permKappaVec(IntegerVector x,
                           IntegerVector y,
                           String method = "unweighted",
                           String alternative = "two.sided",
                           int n_iter = 1000) {

  /// results container

  NumericVector result(6, NA_REAL);

  result.names() =
    CharacterVector({"n", "kappa",
                    "iter_number", "h0_number", "h1_number",
                    "p_value"});

  // kappa for the genuine vector pair

  NumericVector kappa_stat = kappaCpp(x, y, method);

  if(R_IsNA(kappa_stat[1])) return result;

  // resamples of the Y vector and calculation
  // of the resampled kappa statistics

  IntegerMatrix yResamp = resampleIntVec(y, n_iter, false);
  NumericVector kappa_resamp(n_iter);

  for(int i = 0; i < n_iter; ++i) kappa_resamp[i] = kappaCpp(x, yResamp(i, _), method)[1];

  kappa_resamp = na_omit(kappa_resamp);
  int n_eff_iter = kappa_resamp.length();

  if(n_eff_iter == 0) return result;

  // comparison of the resampled correlation coefficient
  // with the correlation coefficient for the genuine vector pair

  NumericVector h0_vector(n_eff_iter, 1.0);

  if(alternative == "two.sided") {

    NumericVector resampAbs = Abs(kappa_resamp, false);

    h0_vector = ifelse(resampAbs < std::abs(kappa_stat[1]), 0.0, 1.0);

  } else if(alternative == "less") {

    h0_vector = ifelse(kappa_resamp > kappa_stat[1], 0.0, 1.0);

  } else {

    h0_vector = ifelse(kappa_resamp < kappa_stat[1], 0.0, 1.0);

  }

  /// calculation of the numbers of cases the H0 and H1 is met,
  /// p value

  double h0_number = sum(h0_vector);

  double h1_number = n_eff_iter - h0_number;

  double p_value;

  if(h0_number != 0) {

    p_value = h0_number/n_eff_iter;

  } else {

    p_value = (h0_number + 1)/n_eff_iter;

  }

  // output

  result[0] = kappa_stat[0];
  result[1] = kappa_stat[1];
  result[2] = n_eff_iter;
  result[3] = h0_number;
  result[4] = h1_number;
  result[5] = p_value;

  return result;

}

// permutation tests for a pair of integer matrices

// [[Rcpp::export]]

NumericMatrix permKappaMtx(IntegerMatrix x,
                           IntegerMatrix y,
                           String method = "unweighted",
                           String alternative = "two.sided",
                           int n_iter = 1000) {

  // result container

  if((x.ncol() != y.ncol()) | (x.nrow() != y.nrow())) {

    stop("Matrices 'x' and 'y' must have the same dimensions");

  }

  int n = x.ncol();

  NumericMatrix result(n, 6);

  colnames(result) =
    CharacterVector({"n",
                    "kappa",
                    "iter_number",
                    "h0_number",
                    "h1_number",
                    "p_value"});


  // kappa permutation tests

  for(int i = 0; i < n; ++i) {

    result(i, _) =
      permKappaVec(x(_, i), y(_, i), method, alternative, n_iter);

  }

  return result;

}

// bootstrap tests for a pair of vectors

// [[Rcpp::export]]

NumericVector bootKappaVec(IntegerVector x,
                           IntegerVector y,
                           String method = "unweighted",
                           String ci_type = "bca",
                           double conf_level = 0.95,
                           int n_iter = 1000) {

  // result container

  NumericVector result(9, NA_REAL);

  result.names() =
    CharacterVector({"n", "kappa",
                    "boot_mean", "lower_ci", "upper_ci",
                    "iter_number", "h0_number", "h1_number",
                    "p_value"});

  // kappa statistic for the genuine vector pair

  NumericVector kappa_stat = kappaCpp(x, y, method);
  int n = x.length();

  //if(R_IsNA(kappa_stat[1])) return result;

  // resamples of the vector pair and calculation
  // of the resampled kappa statistics coefficients

  NumericVector idx = intSeq(0, n - 1);

  NumericMatrix idx_resamp = resampleVec(as<NumericVector>(idx), n_iter, true);

  NumericVector kappa_resamp(n_iter);

  for(int i = 0; i < n_iter; ++i) {

    kappa_resamp[i] =
      kappaCpp(x[idx_resamp(i, _)], y[idx_resamp(i, _)], method)[1];

  }

  kappa_resamp = na_omit(kappa_resamp);
  int n_eff_iter = kappa_resamp.length();

  if(n_eff_iter == 0) return result;

  // bootstrapped mean and confidence intervals

  double bootMean = mean(kappa_resamp);

  NumericVector CI(2);

  if(ci_type == "bca") {

    CI = BCA(kappa_resamp, conf_level);

  } else {

    CI = perCI(kappa_resamp, conf_level);

  }

  // significance: numbers of events in support of the H0 hypothesis,
  // i.e. that the bootstrapped mean is 0.

  NumericVector h0_vector(n_eff_iter);

  h0_vector = ifelse(kappa_resamp <= 0.0, 1.0, 0.0);

  double h0_number = sum(h0_vector);
  double h1_number = n_eff_iter - h0_number;

  double p_value;

  if(h0_number != 0) {

    p_value = h0_number/n_eff_iter;

  } else {

    p_value = (h0_number + 1)/n_eff_iter;

  }

  // output

  result[0] = kappa_stat[0];
  result[1] = kappa_stat[1];
  result[2] = bootMean;
  result[3] = CI[0];
  result[4] = CI[1];
  result[5] = n_eff_iter;
  result[6] = h0_number;
  result[7] = h1_number;
  result[8] = p_value;

  return result;

}

// bootstrap tests for a pair of matrices

// [[Rcpp::export]]

NumericMatrix bootKappaMtx(IntegerMatrix x,
                           IntegerMatrix y,
                           String method = "unweighted",
                           String ci_type = "bca",
                           double conf_level = 0.95,
                           int n_iter = 1000) {

  // result container

  if((x.ncol() != y.ncol()) | (x.nrow() != y.nrow())) {

    stop("Matrices 'x' and 'y' must have the same dimensions");

  }

  int n = x.ncol();

  NumericMatrix result(n, 9);

  colnames(result) =
    CharacterVector({"n",
                    "kappa",
                    "boot_mean",
                    "lower_ci",
                    "upper_ci",
                    "iter_number",
                    "h0_number",
                    "h1_number",
                    "p_value"});

  // Kappa bootstrap tests

  for(int i = 0; i < n; ++i) {

    result(i, _) =
      bootKappaVec(x(_, i), y(_, i), method, ci_type, conf_level, n_iter);

  }

  return result;

}

// END
