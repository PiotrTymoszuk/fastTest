/*** Two-sample Kolmogorov-Smirnov tests*/

#include <Rcpp.h>
#include <Rmath.h>
#include "vectorUtils.h"
#include "transformUtils.h"
#include "numericUtils.h"

using namespace Rcpp;

// Kolmogorov-Smirnov test for a pair of vectors
// NA's are silently removed. Asymptotic p values

//[[Rcpp::export]]

NumericVector ksTestCpp(NumericVector x,
                        NumericVector y,
                        String alternative = "two.sided") {

  /// entry control and processing of the vectors

  x = na_omit(x);
  y = na_omit(y);

  double n_x = x.length() * 1.0;
  double n_y = y.length() * 1.0;
  double ties = 0.0;

  if(n_x < 2) stop("At least two complete observations in 'x' required");
  if(n_y < 2) stop("At least two complete observations in 'y' required");

  /// calculating the distribution difference and KS test statistic

  NumericVector w = Concat(x, y);

  NumericVector w_ranked = ifelse(order_(w) <= n_x, 1/n_x, -1/n_y);

  NumericVector z = cumsum(w_ranked);

  if((unique(w).length() * 1.0) < (n_x + n_y)) {

    /// if there are ties the testing statistic will be computed
    /// only for positions of Z corresponding to unique values

    ties = 1.0;

    IntegerVector z_idx = seq(0, z.length() - 1);
    LogicalVector z_unique_lgl = ifelse(diff(w.sort()) != 0, true, false);
    z_unique_lgl.push_back(true);

    z = z[z_unique_lgl];

  }

  double d;

  if(alternative == "two.sided") {

    d = max(Abs(z, false));

  } else if(alternative == "greater") {

    d = max(z);

  } else {

    d = -min(z);

  }

  // p value and the output

  double p_value = pSmirnov(d, n_x, n_y, alternative)[0];

  NumericVector res = NumericVector::create(n_x,
                                            n_y,
                                            ties,
                                            d,
                                            p_value);

  res.names() = CharacterVector({"n1", "n2", "ties", "d", "p_value"});

  return res;

}

// Kolmogorov-Smirnov test for a numeric vector with
// an integer splitting factor

//[[Rcpp::export]]

NumericVector ksTestVec(NumericVector x,
                        IntegerVector f,
                        String alternative = "two.sided") {

  List x_splits = Split(x, f);

  if(x_splits.length() < 2) stop("Only one strata available");

  return ksTestCpp(x_splits[0], x_splits[1], alternative);

}

// Kolmogorov-Smirnov test for a numeric matrix and an integer
// splitting factor

//[[Rcpp::export]]

NumericMatrix ksTestMtx(NumericMatrix x,
                        IntegerVector f,
                        String alternative = "two.sided") {

  // testing

  int xColSize = x.ncol();

  NumericMatrix resMtx(xColSize, 5);

  colnames(resMtx) =
    CharacterVector({"n1", "n2", "ties", "d", "p_value"});

  for(int i = 0; i < xColSize; ++i) {

    NumericVector testInput = x(_, i);

    NumericVector testRes =
      ksTestVec(testInput, f, alternative);

    resMtx(i, _) = testRes;

  }

  // naming and output

  LogicalVector nameCheck = checkNames(x);

  if(nameCheck[1]) rownames(resMtx) = colnames(x);

  return resMtx;

}

// ... and for two numeric matrices, column and row names management
// and extended checks are done by the R function

//[[Rcpp::export]]

NumericMatrix ksTest2Mtx(NumericMatrix x,
                         NumericMatrix y,
                         String alternative = "two.sided") {

  // entry control and result container

  if(x.ncol() != y.ncol()) {

    stop("Matrices 'x' and 'y' must have the same column number");

  }

  int n = x.ncol();

  NumericMatrix resMtx(n, 5);

  colnames(resMtx) =
    CharacterVector({"n1", "n2", "ties", "d", "p_value"});

  for(int i = 0; i < n; ++i) {

    resMtx(i, _) = ksTestCpp(x(_, i), y(_, i), alternative);

  }

  return resMtx;

}

// END
