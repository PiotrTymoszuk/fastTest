/*** Covariance and correlation coefficients */

#include <Rcpp.h>
#include <Rmath.h>
#include <random>

#include "numericUtils.h"
#include "vectorUtils.h"
#include "transformUtils.h"
#include "rankUtils.h"
#include "covariance.h"

using namespace Rcpp;

// covariance for a pair of vectors: Pearson and Spearman

// [[Rcpp::export]]

NumericVector Cov(NumericVector x,
                  NumericVector y,
                  String method = "pearson") {

  // result container

  NumericVector result(2, NA_REAL);

  result.names() =
    CharacterVector::create("n", "cov");

  // entry control: its accomplished by a direct code instead of
  // the `procVectors` function for speed

  if(x.length() != y.length()) stop("Incompatible vector lengths");

  LogicalVector indexes = !(is_na(x) | is_na(y));

  x = x[indexes];
  y = y[indexes];

  int n = x.length();

  result[0] = n;

  if(n < 2) return result;

  if(method != "pearson") {

    x = Rank(x, "average", false);
    y = Rank(y, "average", false);

  }

  // covariance coefficient

  double mean_x = mean(x);
  double mean_y = mean(y);
  double cov_sum = 0.0;

  for (int i = 0; i < n; ++i) cov_sum += (x[i] - mean_x) * (y[i] - mean_y);

  result[1] = cov_sum / (n - 1); // Sample covariance

  return result;

}

// covariance for a matrix

// [[Rcpp::export]]

NumericMatrix CovMtx(NumericMatrix x,
                     String method = "pearson") {

  // result container

  int n = x.ncol();

  NumericMatrix result(n, n);

  if(checkNames(x)[1]) {

    rownames(result) = colnames(x);
    colnames(result) = colnames(x);

  }

  // covariance

  for(int i = 0; i < n; ++i) {

    for(int j = 0; j < n; ++j) {

      result(i, j) = Cov(x(_, i), x(_, j), method)[1];

    }

  }

  return result;

}

// Pearson's and Spearman's correlation coefficients

double rRho(NumericVector x, NumericVector y) {

  // result container

  double result = NA_REAL;

  // Covariance

  NumericVector covar = Cov(x, y, "pearson");

  if(NumericVector::is_na(covar[1])) return result;

  double x_sd = SD(x, true);
  double y_sd = SD(y, true);

  if((x_sd == 0.0) | (y_sd == 0.0)) return result;

  result = covar[1]/(x_sd * y_sd);

  return result;

}

// concordant and discordant pairs for Kandall's tau coefficients

NumericVector countPairs(NumericVector x, NumericVector y) {

  // variable containers
  // P: stores counts of concordant pairs
  // Q: stores counts of discordant pairs
  // X0 and Y0: store numbers of ties

  int n = x.length();

  double P = 0.0, Q = 0.0, X0 = 0.0, Y0 = 0.0;

  for (int i = 0; i < n; ++i) {

    for (int j = i + 1; j < n; ++j) {

      if ((x[i] < x[j]) & (y[i] < y[j])) {

        P += 1.0;

      } else if ((x[i] > x[j]) & (y[i] > y[j])) {

        P += 1.0;

      } else if ((x[i] != x[j]) & (y[i] != y[j])) {

        Q += 1.0;

      } else {

        if (x[i] == x[j]) Y0 += 1.0;
        if (y[i] == y[j]) X0 += 1.0;

      }

    }

  }

  NumericVector result =  NumericVector::create(P, Q, X0, Y0);

  result.names() = CharacterVector::create("P", "Q", "X0", "Y0");

  return result;

}

// Kendall's tauA and tauB correlation coefficients

double tauA(NumericVector x, NumericVector y) {

  NumericVector pairs = countPairs(x, y);

  int n = x.length();

  double n0 = 1.0 * n * (n - 1)/2; // total N of pairs

  return (pairs[0] - pairs[1])/n0;

}

double tauB(NumericVector x, NumericVector y) {

  NumericVector pairs = countPairs(x, y);

  int n = x.length();

  double n0 = 1.0 * n * (n - 1)/2; // total N of pairs

  double denom = (n0 - pairs[2]) * (n0 - pairs[3]);

  if(denom == 0.0) return NA_REAL;

  return (pairs[0] - pairs[1])/std::sqrt(denom);

}

// Xi family of correlation coefficients: xiA calculates
// the correlation coefficients without tie assumption (like TauA),
// xiB calculates the coefficient assuming a possibility that ties exist

double xiA(NumericVector x, NumericVector y) {

  int n = x.length();

  y = orderVector(y, x);

  NumericVector r = Rank(y, "random", false);

  double deltaSum = 0.0;

  for(int i = 0; i < n - 1; ++i) deltaSum += std::abs(r[i + 1] - r[i]);

  double xi = 1 - 3 * deltaSum/(std::pow(n, 2.0) + 1);

  return xi;

}

double xiB(NumericVector x, NumericVector y) {

  int n = x.length();

  y = orderVector(y, x);

  NumericVector r = Rank(y, "random", false);
  NumericVector l = Rank(y, "max", false);

  double deltaSum = 0.0;

  for(int i = 0; i < n - 1; ++i) deltaSum += std::abs(r[i + 1] - r[i]);

  double lSum = 0.0;

  for(int i = 0; i < n; ++i) lSum += l[i] * (n - l[i]);

  double xi = 1 - n * deltaSum/(2 * lSum);

  return xi;

}

// Correlation coefficients for a pair of vectors

// [[Rcpp::export]]

NumericVector Cor(NumericVector x,
                  NumericVector y,
                  String method = "pearson") {

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

  NumericVector result(2, NA_REAL);

  result.names() =
    CharacterVector::create("n", coef_name);

  // entry control: manual code for speed

  if(x.length() != y.length()) stop("Incompatible vector lengths");

  LogicalVector indexes = !(is_na(x) | is_na(y));

  x = x[indexes];
  y = y[indexes];

  int n = x.length();

  result[0] = n;

  if(n < 3) return result;

  if(method == "spearman") {

    x = Rank(x, "average", false);
    y = Rank(y, "average", false);

  }

  // correlation coefficient

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

  result[1] = corFun(x, y);

  return result;

}

// Correlation coefficients for a matrix

// [[Rcpp::export]]

NumericMatrix CorMtx(NumericMatrix x,
                     String method = "pearson") {

  // result container

  int n = x.ncol();

  NumericMatrix result(n, n);

  if(checkNames(x)[1]) {

    rownames(result) = colnames(x);
    colnames(result) = colnames(x);

  }

  // correlation coefficients

  for(int i = 0; i < n; ++i) {

    for(int j = 0; j < n; ++j) {

      result(i, j) = Cor(x(_, i), x(_, j), method)[1];

    }

  }

  return result;

}

// END
