/*** Resampling correlation tests */

#include <Rcpp.h>
#include <Rmath.h>
#include <random>

#include "numericUtils.h"
#include "vectorUtils.h"
#include "transformUtils.h"
#include "rankUtils.h"
#include "covariance.h"

using namespace Rcpp;

// permutation tests for a pair of vectors

// [[Rcpp::export]]

NumericVector permCorVec(NumericVector x,
                         NumericVector y,
                         String method = "pearson",
                         String alternative = "two.sided",
                         int n_iter = 1000) {

  // result container

  String coef_name = "r";

  if(method == "spearman") {

    coef_name = "rho";

  } else if(method == "kendallA") {

    coef_name = "tauA";

  } else if(method == "kendallB") {

    coef_name = "tauB";

  } else if(method == "xiA") {

    coef_name = "xiA";

  } else if(method == "xiB") {

    coef_name = "xiB";

  }

  NumericVector result(6, NA_REAL);

  result.names() =
    CharacterVector::create("n", coef_name,
                            "iter_number", "h0_number", "h1_number",
                            "p_value");

  // entry control: manual code for speed

  if(x.length() != y.length()) stop("Incompatible vector lengths");

  LogicalVector indexes = !(is_na(x) | is_na(y));

  x = x[indexes];
  y = y[indexes];

  int n = x.length();

  if(n < 3) return result;

  if(method == "spearman") {

    x = Rank(x, "average", false);
    y = Rank(y, "average", false);

  }

  // correlation-calculating functions

  double (*corFun)(NumericVector, NumericVector);

  if((method == "pearson") | (method == "spearman")) {

    corFun = rRho;

  } else if(method == "kendallA") {

    corFun = tauA;

  } else if(method == "kendallB") {

    corFun = tauB;

  } else if(method == "xiA") {

    corFun = xiA;

  } else {

    corFun = xiB;

  }

  // correlation coefficient for the genuine vector pair

  double corCoef = corFun(x, y);

  if(R_IsNA(corCoef)) return result;

  // resamples of the Y vector and calculation
  // of the resampled correlation coefficients

  NumericMatrix yResamp = resampleVec(y, n_iter, false);
  NumericVector corResamp(n_iter);

  for(int i = 0; i < n_iter; ++i) corResamp[i] = corFun(x, yResamp(i, _));

  corResamp = na_omit(corResamp);
  int n_eff_iter = corResamp.length();

  if(n_eff_iter == 0) return result;

  // comparison of the resampled correlation coefficient
  // with the correlation coefficient for the genuine vector pair

  NumericVector h0_vector(n_eff_iter, 1.0);

  if(alternative == "two.sided") {

    NumericVector resampAbs = Abs(corResamp, false);

    h0_vector = ifelse(resampAbs < std::abs(corCoef), 0.0, 1.0);

  } else if(alternative == "less") {

    h0_vector = ifelse(corResamp > corCoef, 0.0, 1.0);

  } else {

    h0_vector = ifelse(corResamp < corCoef, 0.0, 1.0);

  }

  // calculation of the numbers of cases the H0 and H1 is met,
  // p value

  double h0_number = sum(h0_vector);

  double h1_number = n_eff_iter - h0_number;

  double p_value;

  if(h0_number != 0) {

    p_value = h0_number/n_eff_iter;

  } else {

    p_value = (h0_number + 1)/n_eff_iter;

  }

  // output

  result[0] = n;
  result[1] = corCoef;
  result[2] = n_eff_iter;
  result[3] = h0_number;
  result[4] = h1_number;
  result[5] = p_value;

  return result;

}

// permutation tests for a matrix

// [[Rcpp::export]]

NumericMatrix permCorMtx(NumericMatrix x,
                         String method = "pearson",
                         String alternative = "two.sided",
                         int n_iter = 1000) {

  // result container

  String coef_name = "r";

  if(method == "spearman") {

    coef_name = "rho";

  } else if(method == "kendallA") {

    coef_name = "tauA";

  } else if(method == "kendallB") {

    coef_name = "tauB";

  } else if(method == "xiA") {

    coef_name = "xiA";

  } else if(method == "xiB") {

    coef_name = "xiB";

  }

  int n = x.ncol();

  NumericMatrix result(n * n, 8);

  colnames(result) =
    CharacterVector::create("variable1", "variable2",
                            "n", coef_name,
                            "iter_number", "h0_number", "h1_number",
                            "p_value");


  // correlation tests

  NumericVector resampleResult(6);

  int resId = 0;

  for(int i = 0; i < n; ++i) {

    for(int j = 0; j < n; ++j) {

      // indexes are coded according to the R scheme

      resampleResult =
        permCorVec(x(_, i), x(_, j), method, alternative, n_iter);

      result(resId, 0) = 1.0 * i + 1.0;
      result(resId, 1) = 1.0 * j + 1.0;

      // filling the result matrix

      for(int k = 0; k < resampleResult.length(); ++k) result(resId, k + 2) = resampleResult[k];

      resId += 1;

    }

  }

  return result;

}

// bootstrap tests for a pair of vectors

// [[Rcpp::export]]

NumericVector bootCorVec(NumericVector x,
                         NumericVector y,
                         String method = "pearson",
                         String ci_type = "bca",
                         double conf_level = 0.95,
                         int n_iter = 1000) {

  // result container

  String coef_name = "r";

  if(method == "spearman") {

    coef_name = "rho";

  } else if(method == "kendallA") {

    coef_name = "tauA";

  } else if(method == "kendallB") {

    coef_name = "tauB";

  } else if(method == "xiA") {

    coef_name = "xiA";

  } else if(method == "xiB") {

    coef_name = "xiB";

  }

  NumericVector result(9, NA_REAL);

  result.names() =
    CharacterVector::create("n", coef_name,
                            "boot_mean", "lower_ci", "upper_ci",
                            "iter_number", "h0_number", "h1_number",
                            "p_value");

  // entry control: manual code for speed

  if(x.length() != y.length()) stop("Incompatible vector lengths");

  LogicalVector indexes = !(is_na(x) | is_na(y));

  x = x[indexes];
  y = y[indexes];

  int n = x.length();

  if(n < 3) return result;

  if(method == "spearman") {

    x = Rank(x, "average", false);
    y = Rank(y, "average", false);

  }

  // correlation-calculating functions

  double (*corFun)(NumericVector, NumericVector);

  if((method == "pearson") | (method == "spearman")) {

    corFun = rRho;

  } else if(method == "kendallA") {

    corFun = tauA;

  } else if(method == "kendallB") {

    corFun = tauB;

  } else if(method == "xiA") {

    corFun = xiA;

  } else {

    corFun = xiB;

  }

  // correlation coefficient for the genuine vector pair

  double corCoef = corFun(x, y);

  if(R_IsNA(corCoef)) return result;

  // resamples of the vector pair and calculation
  // of the resampled correlation coefficients

  IntegerVector idx = seq(0, n - 1);

  NumericMatrix idxResamp = resampleVec(as<NumericVector>(idx), n_iter, true);

  NumericVector corResamp(n_iter);

  for(int i = 0; i < n_iter; ++i) {

    corResamp[i] = corFun(x[idxResamp(i, _)], y[idxResamp(i, _)]);

  }

  corResamp = na_omit(corResamp);
  int n_eff_iter = corResamp.length();

  if(n_eff_iter == 0) return result;

  // bootstrapped mean and confidence intervals

  double bootMean = mean(corResamp);

  NumericVector CI(2);

  if(ci_type == "bca") {

    CI = BCA(corResamp, conf_level);

  } else {

    CI = perCI(corResamp, conf_level);

  }

  // significance: numbers of events in support of the H0 hypothesis,
  // i.e. that the bootstrapped mean is 0.

  NumericVector h0_vector(n_eff_iter);

  if(corCoef > 0.0) {

    h0_vector = ifelse(corResamp <= 0.0, 1.0, 0.0);

  } else {

    h0_vector = ifelse(corResamp >= 0.0, 1.0, 0.0);

  }

  double h0_number = sum(h0_vector);
  double h1_number = n_eff_iter - h0_number;

  double p_value;

  if(h0_number != 0) {

    p_value = h0_number/n_eff_iter;

  } else {

    p_value = (h0_number + 1)/n_eff_iter;

  }

  // output

  result[0] = n;
  result[1] = corCoef;
  result[2] = bootMean;
  result[3] = CI[0];
  result[4] = CI[1];
  result[5] = n_eff_iter;
  result[6] = h0_number;
  result[7] = h1_number;
  result[8] = p_value;

  return result;

}

// bootstrap tests for a numeric matrix

// [[Rcpp::export]]

NumericVector bootCorMtx(NumericMatrix x,
                         String method = "pearson",
                         String ci_type = "bca",
                         double conf_level = 0.95,
                         int n_iter = 1000) {

  // result container

  String coef_name = "r";

  if(method == "spearman") {

    coef_name = "rho";

  } else if(method == "kendallA") {

    coef_name = "tauA";

  } else if(method == "kendallB") {

    coef_name = "tauB";

  } else if(method == "xiA") {

    coef_name = "xiA";

  } else if(method == "xiB") {

    coef_name = "xiB";

  }

  int n = x.ncol();

  NumericMatrix result(n * n, 11);

  colnames(result) =
    CharacterVector::create("variable1", "variable2",
                            "n", coef_name,
                            "boot_mean", "lower_ci", "upper_ci",
                            "iter_number", "h0_number", "h1_number",
                            "p_value");

  // correlation tests

  NumericVector resampleResult(9);

  int resId = 0;

  for(int i = 0; i < n; ++i) {

    for(int j = 0; j < n; ++j) {

      // indexes of the tested variables are coded according to the R scheme

      resampleResult =
        bootCorVec(x(_, i), x(_, j), method, ci_type, conf_level, n_iter);

      result(resId, 0) = 1.0 * i + 1.0;
      result(resId, 1) = 1.0 * j + 1.0;

      // filling the result matrix

      for(int k = 0; k < resampleResult.length(); ++k) result(resId, k + 2) = resampleResult[k];

      resId += 1;

    }

  }

  return result;

}

// END
