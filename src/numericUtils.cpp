/*** Numeric statistics */

#include <Rcpp.h>
#include <Rmath.h>
#include "numericUtils.h"
#include "vectorUtils.h"

using namespace Rcpp;

// calculation of median, variance, SD, and SEM

// [[Rcpp::export]]

double Median(NumericVector x, bool na_rm = true) {

  if(na_rm) x = na_omit(x);

  int x_size = x.size();
  double med;

  x.sort();

  if(x_size % 2 != 0) {

    med = x[x_size/2];

  } else {

    med = (x[x_size/2 - 1] + x[x_size/2])/2;

  }

  return med;

}

// [[Rcpp::export]]

double Var(NumericVector x, bool na_rm = true) {

  if(na_rm) x = na_omit(x);

  int n = x.size();
  double x_mean = 0;
  double sum = 0;

  x_mean = mean(x);

  for (int i = 0; i < n; ++i) {

    sum += std::pow(x[i] - x_mean, 2);

  }

  return sum / (n - 1);

}

// [[Rcpp::export]]

double SD(NumericVector x, bool na_rm = true) {

  return std::sqrt(Var(x, na_rm));

}

// [[Rcpp::export]]

double SEM(NumericVector x, bool na_rm = true) {

  return SD(x, na_rm)/std::sqrt(x.length());

}

// calculation of pooled variances of two vectors, which will be used
// by functions computing variants of Cohen's d

double poolVarPaired(NumericVector x, NumericVector y) {

  // as per definition of the paired test, x and y has to be of equal lengths
  // the quality check is done by an upstream function

  int xLen = x.length();

  NumericVector res(xLen);

  for(int i = 0; i < xLen; ++i) {

    res[i] = x[i] - y[i];

  }

  return Var(res);

}

double poolVarWelch(NumericVector x, NumericVector y) {

  x = na_omit(x);
  y = na_omit(y);

  return (Var(x) + Var(y))/2;

}

double poolVarStandard(NumericVector x, NumericVector y) {

  x = na_omit(x);
  y = na_omit(y);

  double mean_x = mean(x);
  double mean_y = mean(y);

  int xLen = x.length();
  int yLen = y.length();

  // denominator of the variance formula

  double n_diff = xLen + yLen - 2;

  if(n_diff <= 0) n_diff = 1;

  // sums of squares

  NumericVector squares_x(xLen);
  NumericVector squares_y(yLen);

  for(int i = 0; i < xLen; ++i) {

    squares_x[i] = x[i] - mean_x;

    squares_x[i] = std::pow(squares_x[i], 2.0);

  }

  for(int j = 0; j < yLen; ++j) {

    squares_y[j] = y[j] - mean_y;

    squares_y[j] = std::pow(squares_y[j], 2.0);

  }

  return (sum(squares_x) + sum(squares_y))/n_diff;

}

// Hodges-Lehman location estimates

// location metrics for an unpaired setting

// [[Rcpp::export]]

double locationStd(NumericVector x, NumericVector y) {

  NumericMatrix outerMtx = outerDelta(x, y, true);

  NumericVector outerVals = matrix2vector(outerMtx, false);

  outerVals = na_omit(outerVals);

  return Median(outerVals, true);

}

// location metrics for a single vector of differences

// [[Rcpp::export]]

double locationPair(NumericVector x) {

  // a modified version of R's code: the matrix of outer sums
  // is square and symmetric, hence the median of non-diagonal elements
  // is supposed to be equal to the median of the upper triangle

  NumericMatrix diffs = outerSum(x, x, true);

  diffs.fill_diag(NA_REAL);

  NumericVector mtxVals = matrix2vector(diffs, false);

  mtxVals = na_omit(mtxVals);

  return Median(mtxVals, true);

}

// computation of quantiles and confidence intervals

// [[Rcpp::export]]

NumericVector Quantile(NumericVector x, NumericVector probs) {

  // calculation of quantiles
  // many thanks to https://github.com/RcppCore/Rcpp/issues/967

  const size_t n = x.size(), np = probs.size();

  if (n == 0) return x;
  if (np == 0) return probs;

  NumericVector index = (n - 1.0) * probs, y = x.sort(), x_hi(np), qs(np);
  NumericVector lo = floor(index), hi = ceiling(index);

  for (size_t i = 0; i < np; ++i) {

    qs[i] = y[lo[i]];
    x_hi[i] = y[hi[i]];

    if ((index[i] > lo[i]) && (x_hi[i] != qs[i])) {

      double h;
      h = index[i] - lo[i];
      qs[i] = (1.- h) * qs[i] + h * x_hi[i];

    }

  }

  return qs;

}

// [[Rcpp::export]]

NumericVector perCI(NumericVector theta, double conf_level = 0.95) {

  // computes percentile confidence intervals

  NumericVector ci_probs{(1 - conf_level)/2, (1 + conf_level)/2};

  return Quantile(theta, ci_probs);

}

// [[Rcpp::export]]

NumericVector BCA(NumericVector theta, double conf_level = 0.95) {

  // computes BCA confidence intervals based on the R code
  // provided by the coxed package
  // https://rdrr.io/cran/coxed/src/R/bca.R

  double low;
  double high;

  low = (1 - conf_level)/2;
  high = 1 - low;

  int sims = theta.size();

  NumericVector low_theta;

  low_theta = ifelse(theta < mean(theta), 1.0, 0.0);

  double z_inv;

  z_inv = sum(low_theta)/sims;

  double z;

  z = R::qnorm(z_inv, 0, 1, 1, 0);

  NumericVector U;

  U = (sims - 1) * (mean(theta) - theta);

  double top;
  double under;
  double a;

  top = sum(pow(U, 3));

  under = sum(pow(U, 2));

  under = 6 * std::pow(under, 1.5);

  a = top/under;

  double lower_inv;
  double upper_inv;

  double q_low;
  double q_high;

  q_low = R::qnorm(low, 0, 1, 1, 0);

  lower_inv = z + (z + q_low)/(1 - a * (z + q_low));

  lower_inv = R::pnorm(lower_inv, 0, 1, 1, 0);

  q_high = R::qnorm(high, 0, 1, 1, 0);

  upper_inv = z + (z + q_high)/(1 - a * (z + q_high));

  upper_inv = R::pnorm(upper_inv, 0, 1, 1, 0);

  NumericVector probs = {lower_inv, upper_inv};

  return Quantile(theta, probs);

}

// END
