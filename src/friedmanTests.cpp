/*** Kruskal-Wallis tests */

#include <Rcpp.h>
#include <Rmath.h>
#include "numericUtils.h"
#include "vectorUtils.h"
#include "transformUtils.h"
#include "rankUtils.h"
#include "contingencyUtils.h"
#include "friedmanTests.h"

using namespace Rcpp;

// a base version for a numeric vector x, an integer vector f defining the
// treatment groups, and an integer vector b, that defines the blocks

// [[Rcpp::export]]

NumericVector friedmanVec(NumericVector x,
                          IntegerVector f,
                          IntegerVector b,
                          bool crash = false) {

  // result container

  NumericVector result(8, NA_REAL);

  result.names() =
    CharacterVector::create("n", "k", "b",
                            "q", "tie_correction",
                            "df", "p_value", "w");

  // entry control and ranking of the vectors

  List procVectors = rankBlockVectors(x, f, b);

  bool redFlag = procVectors[1];

  if(redFlag & crash) stop("Calculation failed");

  if(redFlag) {

    warning("Calculation failed");

    return result;

  }

  NumericMatrix rankMtx = procVectors[0];

  // column-wise sum of ranks and calculation of the test statistic Q

  double k = rankMtx.ncol() * 1.0; // number of treatment groups
  double B = rankMtx.nrow() * 1.0; // number of blocks
  double N = k * B; // number of all observations

  NumericVector colSums(k); // will store sums of ranks for the treatments

  for(int i = 0; i < k; ++i) colSums[i] = sum(rankMtx(_, i));

  NumericVector colSquares = colSums * colSums;

  double Q = 12/(B * k * (k + 1)) * sum(colSquares) - 3 * B * (k + 1);

  // row-wise tie correction

  NumericVector rowTies(B);

  for(int i = 0; i < B; ++i) {

    NumericVector ties = as<NumericVector>(Table(rankMtx(i, _)));

    rowTies[i] = sum(ties * ties * ties - ties);

  }

  double tie_correction = 1 - sum(rowTies)/(B * k * (std::pow(k, 2.0) - 1));

  Q = Q/tie_correction;

  // p-value and effect size W metric

  double df = (k - 1) * 1.0;

  double p_value = R::pchisq(Q, df, false, false);

  if(p_value > 1.0) p_value = 1.0;

  double W = Q/(1.0 * N * (k - 1));

  // output

  result[0] = N;
  result[1] = k;
  result[2] = B;

  result[3] = Q;
  result[4] = tie_correction;
  result[5] = df;

  result[6] = p_value;
  result[7] = W;

  return result;

}

// and a version for the matrix input

// [[Rcpp::export]]

NumericMatrix friedmanMtx(NumericMatrix x,
                          IntegerVector f,
                          IntegerVector b,
                          bool crash = false) {


  int nVars = x.ncol();

  NumericMatrix result(nVars, 8);

  colnames(result) =
    CharacterVector::create("n", "k", "b",
                            "q", "tie_correction",
                            "df", "p_value", "w");

  NumericVector tstResult(8);

  for(int i = 0; i < nVars; ++i) {

    tstResult = friedmanVec(x(_, i), f, b, crash);

    result(i, _) = tstResult;

  }

  if(checkNames(x)[1]) rownames(result) = colnames(x);

  return result;

}

// END
