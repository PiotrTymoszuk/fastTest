/*** Kruskal-Wallis tests */

#include <Rcpp.h>
#include <Rmath.h>
#include "numericUtils.h"
#include "vectorUtils.h"
#include "transformUtils.h"
#include "rankUtils.h"
#include "contingencyUtils.h"
#include "kruskalTests.h"

using namespace Rcpp;

// a base version for a numeric vector x and a splitting factor f
// the tie-correction factor is spacified by:
// https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_test

// [[Rcpp::export]]

NumericVector kruskalVec(NumericVector x,
                         IntegerVector f,
                         bool crash = true) {

  // result container

  NumericVector result(7, NA_REAL);

  result.names() =
    CharacterVector::create("n", "k",
                            "h", "tie_correction",
                            "df", "p_value", "etasq");

  // processing of the vectors and entry control

  List procVectors = rankVectors(x, f);

  bool redFlags = procVectors[2];

  if(redFlags & crash) stop("Calculation failed");

  if(redFlags) {

    warning("Crititcal issues of the input vectors");

    return result;

  }

  NumericVector r = procVectors[0]; // ranks of x
  f = procVectors[1];

  double N = r.length() * 1.0; // the total number of complete observations

  // ties and tie-correction factor

  NumericVector ties = as<NumericVector>(Table(r));

  double tie_correction =
    1 - sum(ties * ties * ties - ties)/(std::pow(N, 3.0) - N);

  // the bare and tie-corrected test statistic H, eta-square effect size

  List r_splits = Split(r, f);
  int k = r_splits.length(); // number of analysis groups

  double groupSquares = 0.0; // will store group mean squares

  for(int i = 0; i < k; ++i) {

    NumericVector group = r_splits[i];

    groupSquares += (1.0 * group.length()) * std::pow(mean(group), 2.0);

  }

  double H = 12/(N * (N + 1)) * groupSquares - 3 * (N + 1);

  H = H/tie_correction;

  double etasq = (H - k + 1)/(N - k);

  if(etasq < 0) etasq = 0.0;

  // p value

  double df = (k - 1) * 1.0;

  double p_value = R::pchisq(H, df, false, false);

  if(p_value > 1.0) p_value = 1.0;

  // output

  result[0] = N;
  result[1] = k;

  result[2] = H;
  result[3] = tie_correction;
  result[4] = df;
  result[5] = p_value;

  result[6] = etasq;

  return result;

}

// a version for a matrix

// [[Rcpp::export]]

NumericVector kruskalMtx(NumericMatrix x,
                         IntegerVector f,
                         bool crash = true) {

  int nVars = x.ncol();

  NumericMatrix result(nVars, 7);

  colnames(result) =
    CharacterVector::create("n", "k",
                            "h", "tie_correction",
                            "df", "p_value", "etasq");

  NumericVector tstResult(7);

  for(int i = 0; i < nVars; ++i) {

    tstResult = kruskalVec(x(_, i), f, crash);

    result(i, _) = tstResult;

  }

  if(checkNames(x)[1]) rownames(result) = colnames(x);

  return result;

}

// END
