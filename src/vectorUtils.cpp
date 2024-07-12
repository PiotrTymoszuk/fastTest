/*** General calculation helpers used by other Cpp functions */

#include <Rcpp.h>
#include <Rmath.h>
#include "numericUtils.h"
#include "vectorUtils.h"
#include "contingencyUtils.h"

using namespace Rcpp;

// difference of two numeric vectors

NumericVector delta(NumericVector x, NumericVector y) {

  int xLen = x.length();
  NumericVector res(xLen);

  if(xLen != y.length()) stop("Incompatible vector lengths");

  for(int i = 0; i < xLen; ++i) {

    res[i] = y[i] - x[i];

  }

  return res;

}

// product of two vectors

NumericVector product(NumericVector x, NumericVector y) {

  int xLen = x.length();
  NumericVector res(xLen);

  if(xLen != y.length()) stop("Incompatible vector lengths");

  for(int i = 0; i < xLen; ++i) {

    res[i] = y[i] * x[i];

  }

  return res;

}

// sign function

NumericVector Sign(NumericVector x, bool na_rm = true) {

  if(na_rm) x = na_omit(x);

  int xLen = x.length();

  NumericVector res(xLen, NA_REAL);

  for(int i = 0; i < xLen; ++i) {

    if(x[i] == 0) {

      res[i] = 0;

    } else {

      res[i] = x[i]/std::abs(x[i]);

    }

  }

  return res;

}

// absolute values of a vector

NumericVector Abs(NumericVector x, bool na_rm = true) {

  if(na_rm) x = na_omit(x);

  int xLen = x.length();

  NumericVector res(xLen, NA_REAL);

  for(int i = 0; i < xLen; ++i) {

    res[i] = std::abs(x[i]);

  }

  return res;

}

// outer sum and difference (following the convention of the package: the
// second argument minus the first one), outer product

// [[Rcpp::export]]

NumericMatrix outerSum(NumericVector x,
                       NumericVector y,
                       bool na_rm = true) {

  if(na_rm) {

    x = na_omit(x);
    y = na_omit(y);

  }

  int xLen = x.length();
  int yLen = y.length();

  NumericMatrix resMtx(xLen, yLen);

  for(int i = 0; i < xLen; ++i) {

    for(int j = 0; j < yLen; ++j) {

      resMtx(i, j) = x[i] + y[j];

    }

  }

  return resMtx;

}

// [[Rcpp::export]]

NumericMatrix outerDelta(NumericVector x,
                         NumericVector y,
                         bool na_rm = true) {

  if(na_rm) {

    x = na_omit(x);
    y = na_omit(y);

  }

  int xLen = x.length();
  int yLen = y.length();

  NumericMatrix resMtx(xLen, yLen);

  for(int i = 0; i < xLen; ++i) {

    for(int j = 0; j < yLen; ++j) {

      resMtx(i, j) = y[j] - x[i];

    }

  }

  return resMtx;

}

// [[Rcpp::export]]

NumericMatrix outerProduct(NumericVector x,
                           NumericVector y,
                           bool na_rm = true) {

  if(na_rm) {

    x = na_omit(x);
    y = na_omit(y);

  }

  int xLen = x.length();
  int yLen = y.length();

  NumericMatrix resMtx(xLen, yLen);

  for(int i = 0; i < xLen; ++i) {

    for(int j = 0; j < yLen; ++j) {

      resMtx(i, j) = x[i] * y[j];

    }

  }

  return resMtx;

}

// values from a diagonal of an entire numeric matrix and diagonal

NumericVector matrix2vector(NumericMatrix x, bool by_row = true) {

  if(by_row) x = transpose(x);

  NumericVector vals(x);

  vals.attr("dim") = R_NilValue;

  return vals;

}

// sorting and ordering

// Order the elements of x by sorting y

NumericVector orderVector(NumericVector x, NumericVector y) {

  // First create a vector of indices
  IntegerVector idx = seq_along(x) - 1;

  // Then sort that vector by the values of y

  std::sort(idx.begin(),
            idx.end(),
            [&](int i, int j){return y[i] < y[j];});

  // And return x in that order

  return x[idx];

}

// a table with duplicate frequency

IntegerVector countMulti(NumericVector x) {

  IntegerVector xTab = Table(x);
  LogicalVector duplIdx = ifelse(xTab> 0, true, false);

  return xTab[duplIdx];

}

// permutation and bootstraps of a vector: resamples are stored as rows

NumericMatrix resampleVec(NumericVector x,
                          int n_iter = 100,
                          bool replace = false) {

  int xLen = x.length();
  NumericVector resample(xLen);
  NumericMatrix result(n_iter, xLen);

  for(int i = 0; i < n_iter; ++i) {

    resample = sample(x, xLen, replace);

    result(i, _) = resample;

  }

  return result;

}

// string sequences

CharacterVector stringSeq(std::string prefix, int start, int end) {

  IntegerVector intSeq = seq(start, end);

  CharacterVector outString(intSeq.begin(), intSeq.end());

  int n = outString.length();

  for(int i = 0; i < n; ++i) prefix + outString[i];

  return outString;

}

// END
